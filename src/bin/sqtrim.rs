use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::parser::ValueSource;
use clap::{Arg, ArgMatches, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::filter::Filter;
use rsearch::label::checked_abundance;
use rsearch::writer::Writer;

/// Options that need quality scores, which cannot be used with FASTA input.
const QUALITY_OPTIONS: [&str; 6] = ["qual", "max-ee", "max-ee-rate", "ascii", "qmin", "qmax"];

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

/// parses a non-negative number of expected errors.
fn non_negative(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(x) if x >= 0.0 => Ok(x),
        _ => Err(format!("`{s}` is not a non-negative number")),
    }
}

/// the trimming the options ask for, with no filter. Nothing is discarded.
fn filter_from(matches: &ArgMatches) -> Filter {
    let length = |name| matches.get_one::<u64>(name).map(|&x| x as usize);
    let offset: u8 = matches
        .get_one::<String>("ascii")
        .map_or(33, |a| a.parse().unwrap());
    Filter {
        strip_left: length("strip-left").unwrap_or(0),
        strip_right: length("strip-right").unwrap_or(0),
        trunc_keep: length("len"),
        trunc_qual: matches.get_one::<u8>("qual").map(|&q| i16::from(q)),
        trunc_ee: matches.get_one::<f64>("max-ee").copied(),
        trunc_rate: matches.get_one::<f64>("max-ee-rate").copied(),
        min_len: 0,
        offset,
        min_score: matches.get_one::<i16>("qmin").copied().unwrap_or(0),
        max_score: matches
            .get_one::<i16>("qmax")
            .copied()
            .unwrap_or(126 - i16::from(offset)),
        ..Filter::default()
    }
}

fn main() {
    let length = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("LEN")
            .value_parser(value_parser!(u64))
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
    let mut command = Command::new("sqtrim")
        .about("Shorten reads: strip their ends and truncate them.")
        .long_about(
            "Shorten reads: strip their ends and truncate them.\n\n\
             The ends are stripped first, then reads are truncated to a length, then \
             base by base at a low score or before the expected errors exceed a limit \
             or a rate. Reads are never discarded, even when emptied; sqfilt filters \
             them. Records are written in the input format, FASTA unwrapped.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTQ or FASTA file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file in the input format, default to stdout"),
        )
        .arg(length(
            "strip-left",
            "Remove this many bases from the 5' end",
        ))
        .arg(length(
            "strip-right",
            "Remove this many bases from the 3' end",
        ))
        .arg(
            length("len", "Truncate to this length; shorter reads are kept")
                .value_parser(value_parser!(u64).range(1..)),
        )
        .arg(
            Arg::new("qual")
                .long("qual")
                .value_name("QUAL")
                .value_parser(value_parser!(u8).range(0..=93))
                .help("Truncate before the first base with a score at or below this"),
        )
        .arg(
            Arg::new("max-ee")
                .long("max-ee")
                .value_name("EE")
                .value_parser(non_negative)
                .help("Truncate before the expected errors exceed this"),
        )
        .arg(
            Arg::new("max-ee-rate")
                .long("max-ee-rate")
                .value_name("RATE")
                .value_parser(non_negative)
                .help("Truncate before the expected errors per base exceed this"),
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

    let filter = filter_from(&matches);
    let offset = i16::from(filter.offset);
    let range_error = if filter.min_score > filter.max_score {
        Some("--qmin cannot be greater than --qmax")
    } else if offset + filter.min_score < 33 {
        Some("the sum of --ascii and --qmin must be at least 33")
    } else if offset + filter.max_score > 126 {
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

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, is_fasta, 0);
    let fail_after = |writer: &mut Writer<Box<dyn Write>>, message: &str| -> ! {
        // the reads before an error are written, as vsearch does
        writer.flush().unwrap();
        fail(message);
    };

    let mut record = Record::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail_after(&mut writer, &e.to_string());
        }
        if record.is_empty() {
            if let Some(warning) = reader.stripped_warning() {
                warn(&warning);
            }
            break;
        }

        let qual = record.is_fastq().then(|| record.qual());
        let outcome = match filter.apply(record.seq(), qual, 1) {
            Ok(outcome) => outcome,
            Err(e) => fail_after(&mut writer, &format!("{e} in read {}", record.id())),
        };
        // vsearch reads, and checks, the abundance of every read, after its
        // scores
        if let Err(e) = checked_abundance(record.header()) {
            fail_after(&mut writer, &format!("{e} in read {}", record.header()));
        }
        let qual = qual.map_or(&[][..], |q| outcome.trim(q));
        writer
            .write(record.header(), None, outcome.trim(record.seq()), qual)
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }

    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
