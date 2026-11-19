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

/// The quality cutoffs of the last table. Reads pass while every score so
/// far is above the cutoff. As the cutoffs ascend, the ones a read still
/// passes are always the first ones.
const Q_CUTOFFS: [i64; 4] = [5, 10, 15, 20];
/// The expected-error cutoffs: reads pass while their expected errors so
/// far are at most the cutoff. Descending, likewise.
const EE_CUTOFFS: [f64; 4] = [1.0, 0.5, 0.25, 0.1];

/// The number of printable symbols, from `!` to `~`.
const PRINTABLE: usize = 94;

/// the position of the first symbol at most `limit`, searched a chunk at a
/// time so that the scan vectorizes.
fn first_at_most(symbols: &[u8], limit: u8) -> Option<usize> {
    const CHUNK: usize = 16;
    let (chunks, rest) = symbols.as_chunks::<CHUNK>();
    for (i, chunk) in chunks.iter().enumerate() {
        if chunk.iter().fold(u8::MAX, |low, &q| low.min(q)) <= limit {
            return chunk
                .iter()
                .position(|&q| q <= limit)
                .map(|at| i * CHUNK + at);
        }
    }
    let done = chunks.len() * CHUNK;
    rest.iter().position(|&q| q <= limit).map(|at| done + at)
}

/// The error probability of a quality score, as vsearch computes it.
fn probability(score: f64) -> f64 {
    10f64.powf(-score / 10.0)
}

/// Counts per read length and per position, as vsearch's `--fastq_stats`.
struct Stats {
    /// The score and error probability of each symbol.
    score: [i64; 256],
    error: [f64; 256],
    /// The offset of scores.
    offset: u8,
    /// The lowest and highest scores accepted.
    min_score: i64,
    max_score: i64,
    /// Reads of each length.
    lengths: Vec<u64>,
    /// Per position, the count of each score symbol and the sum over reads
    /// of the expected errors up to it.
    symbols: Vec<[u64; 256]>,
    /// The latest counts of score symbols, which are all printable (from
    /// `!` on), in narrow counters that fit the cache; added to `symbols`
    /// before they can overflow.
    pending: Vec<[u16; PRINTABLE]>,
    pending_reads: u16,
    sum_ee: Vec<f64>,
    /// Per position and cutoff, the reads whose last passing position it
    /// is; turned into the reads passing there by `spread`.
    ee_pass: Vec<[u64; 4]>,
    q_pass: Vec<[u64; 4]>,
}

impl Stats {
    fn new(offset: u8, min_score: i64, max_score: i64) -> Stats {
        let mut score = [0; 256];
        let mut error = [0.0; 256];
        for symbol in 0..256 {
            score[symbol] = symbol as i64 - i64::from(offset);
            error[symbol] = probability(score[symbol] as f64);
        }
        Stats {
            score,
            error,
            offset,
            min_score,
            max_score,
            lengths: Vec::new(),
            symbols: Vec::new(),
            pending: Vec::new(),
            pending_reads: 0,
            sum_ee: Vec::new(),
            ee_pass: Vec::new(),
            q_pass: Vec::new(),
        }
    }

    /// checks the lowest and highest score of a read against the range
    /// accepted, as vsearch does before counting it.
    fn check(&self, qual: &[u8]) -> Result<(), QualityError> {
        if qual.is_empty() {
            return Ok(());
        }
        // a fold over values, which vectorizes where min() and max() over
        // references do not
        let (lowest, highest) = qual
            .iter()
            .fold((u8::MAX, 0), |(low, high), &q| (low.min(q), high.max(q)));
        for symbol in [lowest, highest] {
            let score = self.score[usize::from(symbol)];
            if score < self.min_score || score > self.max_score {
                return Err(QualityError {
                    score: score as i16,
                    min_score: self.min_score as i16,
                    max_score: self.max_score as i16,
                });
            }
        }
        Ok(())
    }

    /// counts a read.
    fn add(&mut self, qual: &[u8]) {
        let length = qual.len();
        if length + 1 > self.lengths.len() {
            self.lengths.resize(length + 1, 0);
            self.symbols.resize(length + 1, [0; 256]);
            self.pending.resize(length + 1, [0; PRINTABLE]);
            self.sum_ee.resize(length + 1, 0.0);
            self.ee_pass.resize(length + 1, [0; 4]);
            self.q_pass.resize(length + 1, [0; 4]);
        }
        self.lengths[length] += 1;

        // three passes, so that the counting, which has no dependencies
        // between positions, is kept apart from the running sums
        if self.pending_reads == u16::MAX {
            self.flush();
        }
        self.pending_reads += 1;
        for (row, &symbol) in self.pending.iter_mut().zip(qual) {
            row[usize::from(symbol - b'!')] += 1;
        }

        // a read passes a quality cutoff up to its first score at or below
        // it. The cutoffs fail in descending order, and each search starts
        // where the one for the cutoff above it stopped
        let mut q_passing = Q_CUTOFFS.len();
        let mut from = 0;
        while q_passing > 0 {
            let cutoff = Q_CUTOFFS[q_passing - 1] + i64::from(self.offset);
            // no symbol is below '!'
            let Ok(limit) = u8::try_from(cutoff.max(0)) else {
                break;
            };
            let Some(at) = first_at_most(&qual[from..], limit) else {
                break;
            };
            from += at;
            q_passing -= 1;
            if from > 0 {
                self.q_pass[from - 1][q_passing] += 1;
            }
        }

        let mut ee = 0.0;
        let mut ee_passing = EE_CUTOFFS.len();
        for (i, (&symbol, sum)) in qual.iter().zip(self.sum_ee.iter_mut()).enumerate() {
            ee += self.error[usize::from(symbol)];
            *sum += ee;
            while ee_passing > 0 && ee > EE_CUTOFFS[ee_passing - 1] {
                ee_passing -= 1;
                if i > 0 {
                    self.ee_pass[i - 1][ee_passing] += 1;
                }
            }
        }
        // the cutoffs still passed at the end pass at every position
        if length > 0 {
            for cutoff in &mut self.q_pass[length - 1][..q_passing] {
                *cutoff += 1;
            }
            for cutoff in &mut self.ee_pass[length - 1][..ee_passing] {
                *cutoff += 1;
            }
        }
    }

    /// adds the pending counts of symbols to the totals.
    fn flush(&mut self) {
        for (total, pending) in self.symbols.iter_mut().zip(&mut self.pending) {
            for (sum, count) in total[usize::from(b'!')..].iter_mut().zip(pending.iter()) {
                *sum += u64::from(*count);
            }
            *pending = [0; PRINTABLE];
        }
        self.pending_reads = 0;
    }

    /// writes the report.
    fn report(mut self) -> String {
        self.flush();
        // a read passing up to a position passes at every one before it
        for table in [&mut self.ee_pass, &mut self.q_pass] {
            let mut passing = [0u64; 4];
            for row in table.iter_mut().rev() {
                for (sum, count) in passing.iter_mut().zip(row.iter()) {
                    *sum += count;
                }
                *row = passing;
            }
        }
        let reads: u64 = self.lengths.iter().sum();
        let n = reads as f64;
        let residues: f64 = self
            .lengths
            .iter()
            .enumerate()
            .fold(0.0, |sum, (length, &count)| {
                sum + (length as u64 * count) as f64
            });
        let min_length = self.lengths.iter().position(|&c| c > 0).unwrap_or(0);
        let max_length = self.lengths.iter().rposition(|&c| c > 0).unwrap_or(0);
        // reads shorter than each length
        let shorter: Vec<u64> = self
            .lengths
            .iter()
            .scan(0, |sum, &count| {
                *sum += count;
                Some(*sum)
            })
            .collect();
        let at_least = |length: usize| {
            let before = if length == 0 { 0 } else { shorter[length - 1] };
            100.0 * (n - before as f64) / n
        };

        let mut out = String::new();
        out.push_str(
            "\nRead length distribution\n      L           N      Pct   AccPct\n\
             -------  ----------  -------  -------\n",
        );
        if reads > 0 {
            for length in (min_length..=max_length).rev() {
                let count = self.lengths[length];
                if count == 0 {
                    continue;
                }
                let mark = if length == max_length { ">=" } else { "  " };
                writeln!(
                    out,
                    "{mark}{length:>5}  {count:>10}   {:>5.1}%   {:>5.1}%",
                    count as f64 * 100.0 / n,
                    at_least(length)
                )
                .unwrap();
            }
        }

        out.push_str(
            "\nQ score distribution\nASCII    Q       Pe           N      Pct   AccPct\n\
             -----  ---  -------  ----------  -------  -------\n",
        );
        let mut by_symbol = [0u64; 256];
        for position in &self.symbols {
            for (sum, count) in by_symbol.iter_mut().zip(position) {
                *sum += count;
            }
        }
        let mut so_far = 0;
        for symbol in (0..256).rev() {
            let count = by_symbol[symbol];
            if count == 0 {
                continue;
            }
            so_far += count;
            writeln!(
                out,
                "    {}  {:>3}  {:>7.5}  {count:>10}  {:>6.1}%  {:>6.1}%",
                char::from(symbol as u8),
                self.score[symbol],
                self.error[symbol],
                100.0 * count as f64 / residues,
                100.0 * so_far as f64 / residues
            )
            .unwrap();
        }

        out.push_str(
            "\n    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct\n\
             -----  -------  ----  -------  --------  -----  ---------  --------\n",
        );
        for length in 2..=max_length {
            let position = &self.symbols[length - 1];
            let count = position.iter().sum::<u64>() as f64;
            let sum_score: i64 = position
                .iter()
                .zip(&self.score)
                .map(|(&c, &s)| c as i64 * s)
                .sum();
            let sum_error = position
                .iter()
                .zip(&self.error)
                .fold(0.0, |sum, (&c, &p)| sum + c as f64 * p);
            let avg_q = sum_score as f64 / count;
            let avg_p = sum_error / count;
            let avg_ee = self.sum_ee[length - 1] / count;
            let rate = avg_ee / length as f64;
            writeln!(
                out,
                "{length:>5}  {:>6.1}%  {avg_q:>4.1}  {:>7.5}  {avg_p:>8.6}  {avg_ee:>5.2}  \
                 {rate:>9.6}  {:>7.3}%",
                at_least(length),
                probability(avg_q),
                100.0 * rate
            )
            .unwrap();
        }

        out.push_str(
            "\n    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000\n\
             -----  -------  -------  -------  -------  -------  -------  -------  -------\n",
        );
        // from the last length where a read still passes the loosest cutoff
        let last = self
            .ee_pass
            .iter()
            .position(|row| row[0] == 0)
            .unwrap_or(self.ee_pass.len());
        for length in (1..=last).rev() {
            let row = self.ee_pass[length - 1];
            write!(out, "{length:>5}").unwrap();
            for count in row {
                write!(out, "  {count:>7}").unwrap();
            }
            for count in row {
                write!(out, "  {:>6.2}%", 100.0 * count as f64 / n).unwrap();
            }
            out.push('\n');
        }

        out.push_str(
            "\nTruncate at first Q\n  Len     Q=5    Q=10    Q=15    Q=20\n\
             -----  ------  ------  ------  ------\n",
        );
        let half = (max_length / 2).max(1);
        for length in (half..=max_length).rev() {
            write!(out, "{length:>5}").unwrap();
            for count in self.q_pass[length - 1] {
                write!(out, "  {:>5.1}%", 100.0 * count as f64 / n).unwrap();
            }
            out.push('\n');
        }

        writeln!(out, "\n{reads:>10}  Recs ({:.1}M), 0 too long", n / 1e6).unwrap();
        if reads > 0 {
            writeln!(out, "{:>10.1}  Avg length", residues / n).unwrap();
        }
        writeln!(out, "{:>9.1}M  Bases", residues / 1e6).unwrap();
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
    let mut command = Command::new("fqstat")
        .about("Report the lengths and quality scores of FASTQ reads.")
        .long_about(
            "Report the lengths and quality scores of FASTQ reads.\n\n\
             The report gives the distribution of read lengths and of quality scores, \
             the average quality and expected errors at each length, and how many reads \
             would pass four expected-error and four quality cutoffs if truncated at \
             each length.",
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
                .help("Output file of the report, default to stdout"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .help("ASCII offset of quality scores, default to 33"),
        )
        .arg(score(
            "qmin",
            "Minimum quality score accepted in the input, default to 0",
        ))
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
        if let Err(e) = stats.check(record.qual()) {
            fail(&format!("{e} in read {}", record.id()));
        }
        stats.add(record.qual());
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
