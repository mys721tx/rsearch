use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fmt::Write as _;
use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::merge::{Merger, Reason, Settings};
use rsearch::writer::Writer;

/// The pairs vsearch reads before merging them.
const CHUNK: usize = 500;

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

fn open(path: &str) -> Box<dyn Read> {
    Box::new(File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}"))))
}

fn create(path: &str) -> Box<dyn Write> {
    Box::new(BufWriter::new(
        File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))),
    ))
}

/// parses a percentage from 0 to 100.
fn percentage(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(value) if (0.0..=100.0).contains(&value) => Ok(value),
        _ => Err("expected a number from 0 to 100".to_string()),
    }
}

/// An expected error, with more decimals the smaller it is, as vsearch's
/// `--eetabbedout`.
fn ee_value(out: &mut String, value: f64) {
    let precision = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
        .iter()
        .position(|&limit| value < limit)
        .map_or(4, |i| 13 - i);
    write!(out, "{value:.precision$}").unwrap();
}

/// The counts and sums of vsearch's merging report.
#[derive(Default)]
struct Report {
    pairs: u64,
    merged: u64,
    failed: [u64; Reason::ALL.len()],
    read_lengths: f64,
    fragment_lengths: f64,
    squared_fragment_lengths: f64,
    ee_forward: f64,
    ee_reverse: f64,
    ee_merged: f64,
    forward_errors: u64,
    reverse_errors: u64,
}

impl Report {
    fn text(&self) -> String {
        let mut out = String::new();
        let percent = |count: u64| {
            if self.pairs > 0 {
                format!(" ({:.1}%)", 100.0 * count as f64 / self.pairs as f64)
            } else {
                String::new()
            }
        };
        let not_merged = self.pairs - self.merged;
        writeln!(out, "{:>10}  Pairs", self.pairs).unwrap();
        writeln!(out, "{:>10}  Merged{}", self.merged, percent(self.merged)).unwrap();
        writeln!(out, "{:>10}  Not merged{}", not_merged, percent(not_merged)).unwrap();
        if not_merged > 0 {
            out.push_str("\nPairs that failed merging due to various reasons:\n");
        }
        for (reason, &count) in Reason::ALL.iter().zip(&self.failed) {
            if count > 0 {
                writeln!(out, "{count:>10}  {}", reason.description()).unwrap();
            }
        }
        out.push('\n');
        if self.pairs > 0 {
            out.push_str("Statistics of all reads:\n");
            let mean = self.read_lengths / (2.0 * self.pairs as f64);
            writeln!(out, "{mean:>10.2}  Mean read length").unwrap();
        }
        if self.merged > 0 {
            let merged = self.merged as f64;
            let mean = self.fragment_lengths / merged;
            let stdev = ((self.squared_fragment_lengths - (2.0 * mean * self.fragment_lengths)
                + (mean * mean * merged))
                / (merged + 0.0))
                .sqrt();
            out.push_str("\nStatistics of merged reads:\n");
            for (value, what) in [
                (mean, "Mean fragment length"),
                (stdev, "Standard deviation of fragment length"),
                (
                    self.ee_forward / merged,
                    "Mean expected error in forward sequences",
                ),
                (
                    self.ee_reverse / merged,
                    "Mean expected error in reverse sequences",
                ),
                (
                    self.ee_merged / merged,
                    "Mean expected error in merged sequences",
                ),
                (
                    self.forward_errors as f64 / merged,
                    "Mean observed errors in merged region of forward sequences",
                ),
                (
                    self.reverse_errors as f64 / merged,
                    "Mean observed errors in merged region of reverse sequences",
                ),
                (
                    (self.forward_errors + self.reverse_errors) as f64 / merged,
                    "Mean observed errors in merged region",
                ),
            ] {
                writeln!(out, "{value:>10.2}  {what}").unwrap();
            }
        }
        out
    }
}

fn main() {
    let file = |name: &'static str, help: &'static str| {
        Arg::new(name).long(name).value_name("FILE").help(help)
    };
    let score = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("QUAL")
            .allow_negative_numbers(true)
            .value_parser(value_parser!(i16))
            .help(help)
    };
    let mut command = Command::new("mergpa")
        .about("Merge overlapping read pairs.")
        .long_about(
            "Merge overlapping read pairs.\n\n\
             Each forward read is merged with the reverse complement of the read at the \
             same position in --reverse where they overlap, found by shared 5-mers and \
             scored by ungapped alignment; a pair merges when one overlap alone scores \
             high enough. Merged bases get posterior scores, and the forward header is \
             kept whole. Pairs not merged go to --not-merged and --rev-not-merged, \
             upper-cased, N scored 0. Both files are FASTQ, with as many reads each.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTQ file of forward reads, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output FASTQ file of merged pairs, default to stdout"),
        )
        .arg(file("reverse", "Input FASTQ file of reverse reads").required(true))
        .arg(file(
            "not-merged",
            "Output file of the forward reads of pairs not merged",
        ))
        .arg(file(
            "rev-not-merged",
            "Output file of the reverse reads of pairs not merged",
        ))
        .arg(file(
            "report",
            "Output file of the counts of pairs merged and not, and why",
        ))
        .arg(file(
            "ee-table",
            "Output file of the expected errors of the overlap in each read and the \
             differences of each merged pair",
        ))
        .arg(
            Arg::new("min-overlap")
                .long("min-overlap")
                .value_name("LEN")
                .value_parser(value_parser!(u32).range(5..))
                .help("Shortest overlap, at least 5, default to 10"),
        )
        .arg(
            Arg::new("max-diffs")
                .long("max-diffs")
                .value_name("COUNT")
                .value_parser(value_parser!(u32))
                .help("Most differences in the overlap, default to 10"),
        )
        .arg(
            Arg::new("max-diff-pct")
                .long("max-diff-pct")
                .value_name("PCT")
                .allow_negative_numbers(true)
                .value_parser(percentage)
                .help("Highest percentage of differences in the overlap, default to 100"),
        )
        .arg(
            Arg::new("stagger")
                .long("stagger")
                .action(ArgAction::SetTrue)
                .help("Merge pairs whose reverse read runs past the start of the forward read"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .help("ASCII offset of quality scores, read and written, default to 33"),
        )
        .arg(score(
            "qmin",
            "Minimum quality score accepted in the input, default to 0",
        ))
        .arg(score(
            "qmax",
            "Maximum quality score accepted in the input, default to 126 - ascii",
        ))
        .arg(score(
            "qmin-out",
            "Minimum quality score of merged bases, default to 0",
        ))
        .arg(score(
            "qmax-out",
            "Maximum quality score of merged bases, default to 41",
        ));
    let matches = command.get_matches_mut();

    let offset: u8 = matches
        .get_one::<String>("ascii")
        .map_or(33, |a| a.parse().unwrap());
    let get = |name: &str, default: i16| matches.get_one::<i16>(name).copied().unwrap_or(default);
    let settings = Settings {
        min_overlap: matches
            .get_one::<u32>("min-overlap")
            .map_or(10, |&n| n as usize),
        max_diffs: matches
            .get_one::<u32>("max-diffs")
            .map_or(10, |&n| n as usize),
        max_diff_pct: matches
            .get_one::<f64>("max-diff-pct")
            .copied()
            .unwrap_or(100.0),
        stagger: matches.get_flag("stagger"),
        offset,
        min_score: get("qmin", 0),
        max_score: get("qmax", 126 - i16::from(offset)),
        min_score_out: get("qmin-out", 0),
        max_score_out: get("qmax-out", 41),
    };
    let ascii = i16::from(offset);
    let range_error = if settings.min_score > settings.max_score {
        Some("--qmin cannot be greater than --qmax")
    } else if ascii + settings.min_score < 33 {
        Some("the sum of --ascii and --qmin must be at least 33")
    } else if ascii + settings.max_score > 126 {
        Some("the sum of --ascii and --qmax must be at most 126")
    } else if settings.min_score_out > settings.max_score_out {
        Some("--qmin-out cannot be greater than --qmax-out")
    } else if ascii + settings.min_score_out < 33 {
        Some("the sum of --ascii and --qmin-out must be at least 33")
    } else if ascii + settings.max_score_out > 126 {
        Some("the sum of --ascii and --qmax-out must be at most 126")
    } else {
        None
    };
    if let Some(message) = range_error {
        command
            .error(clap::error::ErrorKind::ValueValidation, message)
            .exit();
    }

    let fin = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => open(f),
        None => Box::new(io::stdin()),
    };
    let mut forward = Reader::new(BufReader::new(fin));
    let mut reverse = Reader::new(BufReader::new(open(
        matches.get_one::<String>("reverse").unwrap(),
    )));
    for reader in [&mut forward, &mut reverse] {
        if reader.format().unwrap_or_else(|e| fail(&e.to_string())) == Some(Format::Fasta) {
            fail("FASTQ file expected, FASTA file found");
        }
    }

    let writer = |out: Box<dyn Write>| Writer::new(out, false, 0);
    let mut merged_out = writer(match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(io::stdout()),
    });
    let output = |name: &str| matches.get_one::<String>(name).map(|f| writer(create(f)));
    let mut not_merged = output("not-merged");
    let mut rev_not_merged = output("rev-not-merged");
    let mut report_out = matches.get_one::<String>("report").map(|f| create(f));
    let mut ee_out = matches.get_one::<String>("ee-table").map(|f| create(f));

    let mut merger = Merger::new(settings);
    let mut report = Report::default();
    let mut ee_line = String::new();
    let cannot_write = |e: io::Error| -> ! { fail(&format!("cannot write: {e}")) };
    // vsearch reads pairs in chunks and checks each chunk as it reads it,
    // before merging any of its pairs. An error in reading comes before an
    // error in the scores of an earlier pair of the same chunk
    let mut chunk: Vec<(Record, Record)> =
        (0..CHUNK).map(|_| (Record::new(), Record::new())).collect();
    let mut finished = false;
    while !finished {
        let mut filled = 0;
        while filled < CHUNK {
            let (record, mate) = &mut chunk[filled];
            if let Err(e) = forward.read(record) {
                fail(&e.to_string());
            }
            if record.is_empty() {
                finished = true;
                break;
            }
            if let Err(e) = reverse.read(mate) {
                fail(&e.to_string());
            }
            if mate.is_empty() {
                fail("More forward reads than reverse reads");
            }
            // vsearch reads, and checks, the abundance of both reads
            for header in [record.header(), mate.header()] {
                if let Err(e) = checked_abundance(header) {
                    fail(&format!("{e} in read {header}"));
                }
            }
            filled += 1;
        }
        for (record, mate) in &chunk[..filled] {
            report.pairs += 1;
            report.read_lengths += (record.seq().len() + mate.seq().len()) as f64;
            let outcome = merger
                .merge((record.seq(), record.qual()), (mate.seq(), mate.qual()))
                .unwrap_or_else(|e| fail(&format!("{e} in read {}", record.id())));
            match outcome {
                None => {
                    let m = merger.merged();
                    report.merged += 1;
                    let length = m.seq.len();
                    report.fragment_lengths += length as f64;
                    report.squared_fragment_lengths += (length * length) as f64;
                    report.ee_merged += m.ee_merged;
                    report.ee_forward += m.ee_forward;
                    report.ee_reverse += m.ee_reverse;
                    report.forward_errors += m.forward_errors as u64;
                    report.reverse_errors += m.reverse_errors as u64;
                    merged_out
                        .write(record.header(), None, &m.seq, &m.qual)
                        .unwrap_or_else(|e| cannot_write(e));
                    if let Some(out) = ee_out.as_mut() {
                        ee_line.clear();
                        ee_value(&mut ee_line, m.ee_forward);
                        ee_line.push('\t');
                        ee_value(&mut ee_line, m.ee_reverse);
                        writeln!(ee_line, "\t{}\t{}", m.forward_errors, m.reverse_errors).unwrap();
                        out.write_all(ee_line.as_bytes())
                            .unwrap_or_else(|e| cannot_write(e));
                    }
                }
                Some(reason) => {
                    let index = Reason::ALL.iter().position(|&r| r == reason).unwrap();
                    report.failed[index] += 1;
                    for (out, header, (seq, qual)) in [
                        (not_merged.as_mut(), record.header(), merger.forward()),
                        (rev_not_merged.as_mut(), mate.header(), merger.reverse()),
                    ] {
                        if let Some(out) = out {
                            out.write(header, None, seq, qual)
                                .unwrap_or_else(|e| cannot_write(e));
                        }
                    }
                }
            }
        }
    }
    let mut mate = Record::new();
    if let Err(e) = reverse.read(&mut mate) {
        fail(&e.to_string());
    }
    if !mate.is_empty() {
        fail("More reverse reads than forward reads");
    }

    for writer in [
        Some(&mut merged_out),
        not_merged.as_mut(),
        rev_not_merged.as_mut(),
    ]
    .into_iter()
    .flatten()
    {
        writer.flush().unwrap_or_else(|e| cannot_write(e));
    }
    if let Some(out) = report_out.as_mut() {
        out.write_all(report.text().as_bytes())
            .unwrap_or_else(|e| cannot_write(e));
    }
    for out in [report_out.as_mut(), ee_out.as_mut()].into_iter().flatten() {
        out.flush().unwrap_or_else(|e| cannot_write(e));
    }
    for warning in [reverse.stripped_warning(), forward.stripped_warning()]
        .into_iter()
        .flatten()
    {
        warn(&warning);
    }
}
