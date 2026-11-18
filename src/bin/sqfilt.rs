use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::parser::ValueSource;
use clap::{Arg, ArgMatches, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::filter::{Filter, non_n_pct, upper_pct};
use rsearch::label::checked_abundance;
use rsearch::writer::Writer;

/// Options that need quality scores, which cannot be used with FASTA input.
const QUALITY_OPTIONS: [&str; 6] = ["max-ee", "max-ee-rate", "min-qual", "ascii", "qmin", "qmax"];

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
    Box::new(File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))))
}

/// parses a non-negative number of expected errors.
fn non_negative(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(x) if x >= 0.0 => Ok(x),
        _ => Err(format!("`{s}` is not a non-negative number")),
    }
}

/// parses a positive number of expected errors.
fn positive(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(x) if x > 0.0 => Ok(x),
        _ => Err(format!("`{s}` is not a positive number")),
    }
}

/// parses a percentage.
fn percentage(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(x) if (0.0..=100.0).contains(&x) => Ok(x),
        _ => Err(format!("`{s}` is not a percentage from 0 to 100")),
    }
}

/// The filters: those `Filter` applies, and the shares of unmasked
/// residues.
struct Filters {
    filter: Filter,
    /// Limits on the percentage of upper-case residues.
    upper: (f64, f64),
    /// Limits on the percentage of residues other than N.
    non_n: (f64, f64),
}

impl Filters {
    fn from(matches: &ArgMatches) -> Filters {
        let length = |name| matches.get_one::<u64>(name).map(|&x| x as usize);
        let pct = |name, default| matches.get_one::<f64>(name).copied().unwrap_or(default);
        let offset: u8 = matches
            .get_one::<String>("ascii")
            .map_or(33, |a| a.parse().unwrap());
        Filters {
            filter: Filter {
                min_len: length("min-len").unwrap_or(1),
                max_len: length("max-len"),
                max_ns: length("max-ns"),
                max_ee: matches.get_one::<f64>("max-ee").copied(),
                max_rate: matches.get_one::<f64>("max-ee-rate").copied(),
                min_qual: matches
                    .get_one::<u8>("min-qual")
                    .map_or(0, |&q| i16::from(q)),
                min_size: matches.get_one::<u64>("min-size").copied().unwrap_or(1),
                max_size: matches.get_one::<u64>("max-size").copied(),
                offset,
                min_score: matches.get_one::<i16>("qmin").copied().unwrap_or(0),
                max_score: matches
                    .get_one::<i16>("qmax")
                    .copied()
                    .unwrap_or(126 - i16::from(offset)),
                ..Filter::default()
            },
            upper: (pct("min-upper-pct", 0.0), pct("max-upper-pct", 100.0)),
            non_n: (pct("min-non-n-pct", 0.0), pct("max-non-n-pct", 100.0)),
        }
    }

    /// whether a read fails a filter; fails on invalid abundances and
    /// scores, as vsearch does.
    fn discards(&self, r: &Record) -> Result<bool, String> {
        let header = r.header();
        // vsearch reads, and checks, every abundance
        let size = checked_abundance(header).map_err(|e| format!("{e} in read {header}"))?;
        let qual = r.is_fastq().then(|| r.qual());
        let outcome = self
            .filter
            .apply(r.seq(), qual, size)
            .map_err(|e| format!("{e} in read {}", r.id()))?;
        let outside = |(min, max): (f64, f64), pct: f64| pct < min || pct > max;
        let masked = |limits: (f64, f64)| limits != (0.0, 100.0);
        Ok(outcome.discarded
            || (masked(self.upper) && outside(self.upper, upper_pct(r.seq())))
            || (masked(self.non_n) && outside(self.non_n, non_n_pct(r.seq()))))
    }
}

/// The four outputs: kept and discarded, forward and reverse.
struct Outputs {
    kept: Writer<Box<dyn Write>>,
    discarded: Option<Writer<Box<dyn Write>>>,
    rev_kept: Option<Writer<Box<dyn Write>>>,
    rev_discarded: Option<Writer<Box<dyn Write>>>,
}

impl Outputs {
    fn flush(&mut self) {
        for writer in [
            Some(&mut self.kept),
            self.discarded.as_mut(),
            self.rev_kept.as_mut(),
            self.rev_discarded.as_mut(),
        ]
        .into_iter()
        .flatten()
        {
            writer
                .flush()
                .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
        }
    }

    /// writes what is buffered and fails with `message`, so that the
    /// records before an error are written, as vsearch does.
    fn fail(&mut self, message: &str) -> ! {
        self.flush();
        fail(message);
    }
}

/// writes a record unchanged, if there is an output.
fn write(out: Option<&mut Writer<Box<dyn Write>>>, r: &Record) {
    if let Some(out) = out {
        out.write(r.header(), None, r.seq(), r.qual())
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
}

fn main() {
    let length = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("LEN")
            .value_parser(value_parser!(u64).range(1..))
            .help(help)
    };
    let pct = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("PCT")
            .value_parser(percentage)
            .help(help)
    };
    let score = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("QUAL")
            .allow_negative_numbers(true)
            .value_parser(value_parser!(i16))
            .help(help)
    };
    let file = |name: &'static str, help: &'static str| {
        Arg::new(name).long(name).value_name("FILE").help(help)
    };
    let mut command = Command::new("sqfilt")
        .about("Keep or discard whole reads, or pairs of reads.")
        .long_about(
            "Keep or discard whole reads, or pairs of reads.\n\n\
             A read is discarded when it fails any filter: length, N bases, expected \
             errors and their rate, lowest score, abundance (;size=), or the shares of \
             upper-case and of non-N residues left by masking. With --reverse, a pair \
             is discarded when either read is. Reads are written unchanged, in the \
             input format, FASTA unwrapped.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTQ or FASTA file (forward reads), default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file of kept reads, default to stdout"),
        )
        .arg(file("discarded", "Output file of discarded reads"))
        .arg(file(
            "reverse",
            "Input file of reverse reads, filtered in pairs with the input",
        ))
        .arg(file("rev-out", "Output file of kept reverse reads").requires("reverse"))
        .arg(file("rev-discarded", "Output file of discarded reverse reads").requires("reverse"))
        .arg(length(
            "min-len",
            "Discard reads shorter than this, default to 1",
        ))
        .arg(length("max-len", "Discard reads longer than this"))
        .arg(
            Arg::new("max-ns")
                .long("max-ns")
                .value_name("COUNT")
                .value_parser(value_parser!(u64))
                .help("Discard reads with more N bases than this"),
        )
        .arg(
            Arg::new("max-ee")
                .long("max-ee")
                .value_name("EE")
                .value_parser(positive)
                .help("Discard reads with more expected errors than this"),
        )
        .arg(
            Arg::new("max-ee-rate")
                .long("max-ee-rate")
                .value_name("RATE")
                .value_parser(non_negative)
                .help("Discard reads with more expected errors per base than this"),
        )
        .arg(
            Arg::new("min-qual")
                .long("min-qual")
                .value_name("QUAL")
                .value_parser(value_parser!(u8).range(0..=93))
                .help("Discard reads with a score below this"),
        )
        .arg(
            Arg::new("min-size")
                .long("min-size")
                .value_name("SIZE")
                .value_parser(value_parser!(u64).range(1..))
                .help("Discard reads with an abundance (;size=) below this"),
        )
        .arg(
            Arg::new("max-size")
                .long("max-size")
                .value_name("SIZE")
                .value_parser(value_parser!(u64).range(1..))
                .help("Discard reads with an abundance (;size=) above this"),
        )
        .arg(pct(
            "min-upper-pct",
            "Discard reads with a lower percentage of upper-case residues",
        ))
        .arg(pct(
            "max-upper-pct",
            "Discard reads with a higher percentage of upper-case residues",
        ))
        .arg(pct(
            "min-non-n-pct",
            "Discard reads with a lower percentage of residues other than N",
        ))
        .arg(pct(
            "max-non-n-pct",
            "Discard reads with a higher percentage of residues other than N",
        ))
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

    let filters = Filters::from(&matches);
    let filter = &filters.filter;
    let offset = i16::from(filter.offset);
    let range_error = if filter.min_score > filter.max_score {
        Some("--qmin cannot be greater than --qmax")
    } else if offset + filter.min_score < 33 {
        Some("the sum of --ascii and --qmin must be at least 33")
    } else if offset + filter.max_score > 126 {
        Some("the sum of --ascii and --qmax must be at most 126")
    } else if filters.upper.0 > filters.upper.1 {
        Some("--min-upper-pct cannot be greater than --max-upper-pct")
    } else if filters.non_n.0 > filters.non_n.1 {
        Some("--min-non-n-pct cannot be greater than --max-non-n-pct")
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
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));

    let mut reverse = matches
        .get_one::<String>("reverse")
        .map(|f| Reader::new(BufReader::new(open(f))));
    if let Some(rev) = reverse.as_mut() {
        let rev_format = rev.format().unwrap_or_else(|e| fail(&e.to_string()));
        if format.is_some() && rev_format.is_some() && rev_format != format {
            fail("forward and reverse reads must be in the same format");
        }
    }

    let is_fasta = format == Some(Format::Fasta);
    if is_fasta
        && let Some(option) = QUALITY_OPTIONS
            .iter()
            .find(|&&id| matches.value_source(id) == Some(ValueSource::CommandLine))
    {
        command
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                format!("--{option} needs quality scores, which FASTA input lacks"),
            )
            .exit();
    }

    let writer = |out: Box<dyn Write>| Writer::new(out, is_fasta, 0);
    let optional = |name: &str| matches.get_one::<String>(name).map(|f| writer(create(f)));
    let mut out = Outputs {
        kept: writer(match matches.get_one::<String>("out").map(String::as_str) {
            Some(f) => create(f),
            None => Box::new(io::stdout()),
        }),
        discarded: optional("discarded"),
        rev_kept: optional("rev-out"),
        rev_discarded: optional("rev-discarded"),
    };

    let (mut r, mut r2) = (Record::new(), Record::new());
    loop {
        if let Err(e) = reader.read(&mut r) {
            out.fail(&e.to_string());
        }
        if let Some(rev) = reverse.as_mut() {
            if let Err(e) = rev.read(&mut r2) {
                out.fail(&e.to_string());
            }
            match (r.is_empty(), r2.is_empty()) {
                (true, false) => out.fail("More reverse reads than forward reads"),
                (false, true) => out.fail("More forward reads than reverse reads"),
                _ => {}
            }
        }
        if r.is_empty() {
            break;
        }

        let discards = |r: &Record| filters.discards(r);
        let discarded = match discards(&r) {
            Ok(d) => d,
            Err(e) => out.fail(&e),
        };
        // a pair is discarded if either read is
        let discarded = discarded
            | match reverse.is_some().then(|| discards(&r2)) {
                Some(Ok(d)) => d,
                Some(Err(e)) => out.fail(&e),
                None => false,
            };

        if discarded {
            write(out.discarded.as_mut(), &r);
            if reverse.is_some() {
                write(out.rev_discarded.as_mut(), &r2);
            }
        } else {
            write(Some(&mut out.kept), &r);
            if reverse.is_some() {
                write(out.rev_kept.as_mut(), &r2);
            }
        }
    }

    out.flush();
    for warning in [Some(&reader), reverse.as_ref()]
        .into_iter()
        .flatten()
        .filter_map(Reader::stripped_warning)
    {
        warn(&warning);
    }
}
