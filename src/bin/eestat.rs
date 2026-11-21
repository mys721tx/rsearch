use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fmt::Write as _;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::filter::QualityError;
use rsearch::report::general;

/// prints a warning, which may hold bytes that are not UTF-8.
fn warn(message: &[u8]) {
    let mut stderr = io::stderr().lock();
    let _ = stderr
        .write_all(b"Warning: ")
        .and_then(|()| stderr.write_all(message))
        .and_then(|()| stderr.write_all(b"\n"));
}

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// The scores a symbol can decode to, negative ones counting as 0: from 0
/// to `~` at offset 33.
const SCORES: usize = 94;
/// Expected errors are binned by thousandths for their quartiles.
const RESOLUTION: i64 = 1000;

/// The error probability of a quality score, as vsearch computes it.
fn probability(score: i64) -> f64 {
    10f64.powf(-score as f64 / 10.0)
}

/// The counts of the expected-error bins at one position, kept sorted and
/// sparse. New bins are buffered and merged in once the buffer is as large
/// as the merged counts.
#[derive(Clone, Default)]
struct Bins {
    counts: Vec<(i64, u64)>,
    pending: Vec<i64>,
}

impl Bins {
    fn add(&mut self, bin: i64) {
        self.pending.push(bin);
        if self.pending.len() >= self.counts.len().max(1024) {
            self.merge();
        }
    }

    fn merge(&mut self) {
        if self.pending.is_empty() {
            return;
        }
        self.pending.sort_unstable();
        let mut merged = Vec::with_capacity(self.counts.len() + self.pending.len());
        let mut old = self.counts.iter().copied().peekable();
        for run in self.pending.chunk_by(|a, b| a == b) {
            let bin = run[0];
            while let Some(&(earlier, count)) = old.peek()
                && earlier < bin
            {
                merged.push((earlier, count));
                old.next();
            }
            let mut count = run.len() as u64;
            if let Some(&(same, before)) = old.peek()
                && same == bin
            {
                count += before;
                old.next();
            }
            merged.push((bin, count));
        }
        merged.extend(old);
        self.counts = merged;
        self.pending.clear();
    }

    /// the bins with their counts, in increasing order.
    fn sorted(&mut self) -> &[(i64, u64)] {
        self.merge();
        &self.counts
    }
}

/// The minimum, quartiles, maximum and mean of values seen in increasing
/// order with their counts, as vsearch computes them. A quartile is the
/// first value whose cumulative count reaches that share of the reads.
struct Summary {
    reads: f64,
    seen: f64,
    sum: f64,
    /// minimum, lower quartile, median, upper quartile, maximum; -1 until
    /// seen.
    order: [f64; 5],
}

impl Summary {
    fn new(reads: u64) -> Summary {
        Summary {
            reads: reads as f64,
            seen: 0.0,
            sum: 0.0,
            order: [-1.0; 5],
        }
    }

    fn observe(&mut self, value: f64, count: f64) {
        self.seen += count;
        self.sum += value * count;
        if self.order[0] < 0.0 {
            self.order[0] = value;
        }
        for (k, share) in [(1, 0.25), (2, 0.50), (3, 0.75)] {
            if self.order[k] < 0.0 && self.seen >= share * self.reads {
                self.order[k] = value;
            }
        }
        self.order[4] = value;
    }

    /// minimum, lower quartile, median, mean, upper quartile, maximum,
    /// each mapped by `map` but the mean, given as `mean` or computed.
    fn columns(&self, map: impl Fn(f64) -> f64, mean: Option<f64>) -> [f64; 6] {
        let [min, low, median, high, max] = self.order.map(map);
        [
            min,
            low,
            median,
            mean.unwrap_or(self.sum / self.reads),
            high,
            max,
        ]
    }
}

/// Counts per position, as vsearch's `--fastq_eestats`.
struct Stats {
    /// The score each symbol decodes to, negative ones as 0, whether it is
    /// accepted, and the error probability of its score, at most 1.
    score: [usize; 128],
    accepted: [bool; 128],
    error: [f64; 128],
    offset: i64,
    min_score: i64,
    max_score: i64,
    reads: u64,
    /// Per position: the reads reaching it, the count of each score, the
    /// sum and the bins of the expected errors up to it.
    reaching: Vec<u64>,
    scores: Vec<[u64; SCORES]>,
    sum_ee: Vec<f64>,
    bins: Vec<Bins>,
}

impl Stats {
    fn new(offset: u8, min_score: i64, max_score: i64) -> Stats {
        let offset = i64::from(offset);
        let mut score = [0; 128];
        let mut accepted = [false; 128];
        let mut error = [0.0; 128];
        for symbol in 0..128 {
            let decoded = symbol as i64 - offset;
            score[symbol] = decoded.max(0) as usize;
            accepted[symbol] = (min_score..=max_score).contains(&decoded);
            error[symbol] = probability(decoded).min(1.0);
        }
        Stats {
            score,
            accepted,
            error,
            offset,
            min_score,
            max_score,
            reads: 0,
            reaching: Vec::new(),
            scores: Vec::new(),
            sum_ee: Vec::new(),
            bins: Vec::new(),
        }
    }

    /// counts a read, or fails at its first score out of range.
    fn add(&mut self, qual: &[u8]) -> Result<(), QualityError> {
        if let Some(&symbol) = qual.iter().find(|&&q| !self.accepted[usize::from(q)]) {
            return Err(QualityError {
                score: (i64::from(symbol) - self.offset) as i16,
                min_score: self.min_score as i16,
                max_score: self.max_score as i16,
            });
        }
        self.reads += 1;
        let length = qual.len();
        if length > self.reaching.len() {
            self.reaching.resize(length, 0);
            self.scores.resize(length, [0; SCORES]);
            self.sum_ee.resize(length, 0.0);
            self.bins.resize(length, Bins::default());
        }
        let mut ee = 0.0;
        for (i, &symbol) in qual.iter().enumerate() {
            let symbol = usize::from(symbol);
            self.reaching[i] += 1;
            self.scores[i][self.score[symbol]] += 1;
            ee += self.error[symbol];
            let bin = (RESOLUTION * (i as i64 + 1)).min((RESOLUTION as f64 * ee) as i64);
            self.bins[i].add(bin);
            self.sum_ee[i] += ee;
        }
        Ok(())
    }

    /// writes the table.
    fn report(mut self) -> String {
        let mut out = String::from(
            "Pos\tRecs\tPctRecs\t\
             Min_Q\tLow_Q\tMed_Q\tMean_Q\tHi_Q\tMax_Q\t\
             Min_Pe\tLow_Pe\tMed_Pe\tMean_Pe\tHi_Pe\tMax_Pe\t\
             Min_EE\tLow_EE\tMed_EE\tMean_EE\tHi_EE\tMax_EE\n",
        );
        for i in 0..self.reaching.len() {
            let reads = self.reaching[i];
            let scores = &self.scores[i];
            let mut quality = Summary::new(reads);
            for (score, &count) in scores.iter().enumerate() {
                if count > 0 {
                    quality.observe(score as f64, count as f64);
                }
            }
            // error probabilities in increasing order: scores decreasing
            let mut error = Summary::new(reads);
            for (score, &count) in scores.iter().enumerate().rev() {
                if count > 0 {
                    error.observe(probability(score as i64), count as f64);
                }
            }
            let mut expected = Summary::new(reads);
            for &(bin, count) in self.bins[i].sorted() {
                expected.observe(bin as f64, count as f64);
            }
            let midpoint = |bin: f64| (bin + 0.5) / RESOLUTION as f64;
            let mean_ee = self.sum_ee[i] / reads as f64;

            write!(
                out,
                "{}\t{reads}\t{:.1}",
                i + 1,
                100.0 * reads as f64 / self.reads as f64
            )
            .unwrap();
            for value in quality.columns(|q| q, None) {
                write!(out, "\t{value:.1}").unwrap();
            }
            for value in error.columns(|p| p, None) {
                write!(out, "\t{}", general(value, 2)).unwrap();
            }
            for value in expected.columns(midpoint, Some(mean_ee)) {
                write!(out, "\t{value:.2}").unwrap();
            }
            out.push('\n');
        }
        out
    }
}

fn main() {
    let score = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("QUAL")
            .allow_negative_numbers(true)
            .value_parser(value_parser!(i16))
            .help(help)
    };
    let mut command = Command::new("eestat")
        .about("Report quality scores and expected errors at each position of FASTQ reads.")
        .long_about(
            "Report quality scores and expected errors at each position of FASTQ reads.\n\n\
             For each position the table gives the reads reaching it, then the minimum, \
             lower quartile, median, mean, upper quartile and maximum of the quality \
             scores, of their error probabilities and of the expected errors up to that \
             position. Negative scores count as zero.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTQ file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file of the table, default to stdout"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .default_value("33")
                .help("ASCII offset of quality scores"),
        )
        .arg(score("qmin", "Minimum quality score accepted in the input").default_value("0"))
        .arg(score(
            "qmax",
            "Maximum quality score accepted in the input, default to 126 - ascii",
        ));
    let matches = command.get_matches_mut();

    let offset: u8 = matches
        .get_one::<String>("ascii")
        .map_or(33, |a| a.parse().unwrap());
    let min_score = matches.get_one::<i16>("qmin").copied().unwrap_or(0);
    let max_score = matches
        .get_one::<i16>("qmax")
        .copied()
        .unwrap_or(126 - i16::from(offset));
    let range_error = if min_score > max_score {
        Some("--qmin cannot be greater than --qmax")
    } else if i16::from(offset) + min_score < 33 {
        Some("the sum of --ascii and --qmin must be at least 33")
    } else if i16::from(offset) + max_score > 126 {
        Some("the sum of --ascii and --qmax must be at most 126")
    } else {
        None
    };
    if let Some(message) = range_error {
        command
            .error(clap::error::ErrorKind::ValueValidation, message)
            .exit();
    }

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));
    if format == Some(Format::Fasta) {
        fail("FASTQ file expected, FASTA file found");
    }

    let mut stats = Stats::new(offset, i64::from(min_score), i64::from(max_score));
    let mut record = Record::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        if let Err(e) = stats.add(record.qual()) {
            fail(&format!("{e} in read {}", record.id()));
        }
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }

    let mut fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    fout.write_all(stats.report().as_bytes())
        .and_then(|()| fout.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
