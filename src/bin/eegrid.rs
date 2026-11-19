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

/// parses an expected-error cutoff, which must be a finite number above 0.
fn cutoff(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(value) if value.is_finite() && value > 0.0 => Ok(value),
        _ => Err("expected a number above 0".to_string()),
    }
}

/// The length cutoffs: from `first` by `step`, up to `last`.
struct Lengths {
    first: u64,
    last: u64,
    step: u64,
}

/// Reads kept by length and expected-error cutoffs, as vsearch's
/// `--fastq_eestats2`.
struct Grid {
    lengths: Lengths,
    max_ee: Vec<f64>,
    /// The error probability of each symbol's score, at most 1, and
    /// whether the score is accepted.
    error: [f64; 128],
    accepted: [bool; 128],
    offset: i64,
    min_score: i64,
    max_score: i64,
    reads: u64,
    residues: u64,
    longest: u64,
    /// The length cutoffs up to the longest read so far, and per cutoff
    /// and expected-error cutoff the reads within both.
    rows: usize,
    counts: Vec<u64>,
}

impl Grid {
    fn new(lengths: Lengths, max_ee: Vec<f64>, offset: u8, min_score: i64, max_score: i64) -> Grid {
        let offset = i64::from(offset);
        let mut error = [0.0; 128];
        let mut accepted = [false; 128];
        for symbol in 0..128 {
            let score = symbol as i64 - offset;
            error[symbol] = 10f64.powf(-score as f64 / 10.0).min(1.0);
            accepted[symbol] = (min_score..=max_score).contains(&score);
        }
        Grid {
            lengths,
            max_ee,
            error,
            accepted,
            offset,
            min_score,
            max_score,
            reads: 0,
            residues: 0,
            longest: 0,
            rows: 0,
            counts: Vec::new(),
        }
    }

    /// counts a read, or fails at its first score out of range.
    fn add(&mut self, qual: &[u8]) -> Result<(), QualityError> {
        // the lowest and highest symbols, in a fold that vectorizes; the
        // first symbol out of range is looked for only when there is one
        let (low, high) = qual
            .iter()
            .fold((u8::MAX, 0), |(low, high), &q| (low.min(q), high.max(q)));
        if !qual.is_empty()
            && !(self.accepted[usize::from(low)] && self.accepted[usize::from(high)])
        {
            let symbol = *qual
                .iter()
                .find(|&&q| !self.accepted[usize::from(q)])
                .unwrap();
            return Err(QualityError {
                score: (i64::from(symbol) - self.offset) as i16,
                min_score: self.min_score as i16,
                max_score: self.max_score as i16,
            });
        }
        self.reads += 1;
        let length = qual.len() as u64;
        self.residues += length;
        if length > self.longest {
            self.longest = length;
            // one row at least, even below the first cutoff
            let high = length.min(self.lengths.last);
            let rows = 1 + (high.saturating_sub(self.lengths.first) / self.lengths.step) as usize;
            if rows > self.rows {
                self.rows = rows;
                self.counts.resize(rows * self.max_ee.len(), 0);
            }
        }

        let mut ee = 0.0;
        let (mut row, mut next) = (0, self.lengths.first);
        for (i, &symbol) in qual.iter().enumerate() {
            ee += self.error[usize::from(symbol)];
            if row < self.rows && i as u64 + 1 == next {
                let counts = &mut self.counts[row * self.max_ee.len()..][..self.max_ee.len()];
                for (count, &cutoff) in counts.iter_mut().zip(&self.max_ee) {
                    if ee <= cutoff {
                        *count += 1;
                    }
                }
                row += 1;
                next += self.lengths.step;
            }
        }
        Ok(())
    }

    /// writes the grid.
    fn report(&self) -> String {
        let mut out = format!("{} reads", self.reads);
        if self.reads > 0 {
            write!(
                out,
                ", max len {}, avg {:.1}",
                self.longest,
                self.residues as f64 / self.reads as f64
            )
            .unwrap();
        }
        out.push_str("\n\nLength");
        for cutoff in &self.max_ee {
            write!(out, "         MaxEE {cutoff:.2}").unwrap();
        }
        out.push_str("\n------");
        out.push_str(&"   ----------------".repeat(self.max_ee.len()));
        out.push('\n');
        for row in 0..self.rows {
            write!(
                out,
                "{:>6}",
                self.lengths.first + row as u64 * self.lengths.step
            )
            .unwrap();
            for &count in &self.counts[row * self.max_ee.len()..][..self.max_ee.len()] {
                write!(
                    out,
                    "   {count:>8}({:>5.1}%)",
                    100.0 * count as f64 / self.reads as f64
                )
                .unwrap();
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
    let length = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("LEN")
            .value_parser(value_parser!(u64).range(1..))
            .help(help)
    };
    let mut command = Command::new("eegrid")
        .about("Report the FASTQ reads kept by length and expected-error cutoffs.")
        .long_about(
            "Report the FASTQ reads kept by length and expected-error cutoffs.\n\n\
             For each length cutoff, from --min-len by --len-step up to the longest read \
             or --max-len, the grid gives the reads at least that long whose expected \
             errors up to it are within each of the --max-ee cutoffs, as if truncated \
             there and filtered.",
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
                .help("Output file of the grid, default to stdout"),
        )
        .arg(length("min-len", "First length cutoff, default to 50"))
        .arg(length(
            "max-len",
            "Last length cutoff, default to the longest read",
        ))
        .arg(length(
            "len-step",
            "Step between length cutoffs, default to 50",
        ))
        .arg(
            Arg::new("max-ee")
                .long("max-ee")
                .value_name("EE,...")
                .value_delimiter(',')
                .value_parser(cutoff)
                .help("Expected-error cutoffs, in the order given, default to 0.5,1,2"),
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

    let lengths = Lengths {
        first: matches.get_one::<u64>("min-len").copied().unwrap_or(50),
        last: matches
            .get_one::<u64>("max-len")
            .copied()
            .unwrap_or(u64::MAX),
        step: matches.get_one::<u64>("len-step").copied().unwrap_or(50),
    };
    let max_ee: Vec<f64> = matches
        .get_many::<f64>("max-ee")
        .map_or_else(|| vec![0.5, 1.0, 2.0], |values| values.copied().collect());
    let offset: u8 = matches
        .get_one::<String>("ascii")
        .map_or(33, |a| a.parse().unwrap());
    let min_score = matches.get_one::<i16>("qmin").copied().unwrap_or(0);
    let max_score = matches
        .get_one::<i16>("qmax")
        .copied()
        .unwrap_or(126 - i16::from(offset));
    let range_error = if lengths.first > lengths.last {
        Some("--min-len cannot be greater than --max-len")
    } else if min_score > max_score {
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

    let mut grid = Grid::new(
        lengths,
        max_ee,
        offset,
        i64::from(min_score),
        i64::from(max_score),
    );
    let mut record = Record::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        if let Err(e) = grid.add(record.qual()) {
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
    fout.write_all(grid.report().as_bytes())
        .and_then(|()| fout.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
