use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::{Label, checked_abundance};
use rsearch::random::{Mt19937_64, draw, seed};
use rsearch::writer::Writer;

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

fn create(path: &str) -> Box<dyn Write> {
    Box::new(File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))))
}

/// parses a percentage above 0 and at most 100.
fn percentage(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(x) if x > 0.0 && x <= 100.0 => Ok(x),
        _ => Err(format!("`{s}` is not a percentage above 0 and at most 100")),
    }
}

/// writes the records with reads in `counts`; with `by_size`, each with its
/// count as its abundance.
fn write_counts(
    writer: &mut Writer<Box<dyn Write>>,
    records: &[(Record, u64)],
    counts: &[u64],
    by_size: bool,
) {
    let label = Label {
        sizeout: true,
        ..Label::default()
    };
    for ((r, _), &count) in records.iter().zip(counts) {
        if count == 0 {
            continue;
        }
        let result = if by_size {
            let header = label.format(r.header(), r.seq(), count, 0, None);
            writer.write(&header, None, r.seq(), r.qual())
        } else {
            writer.write(r.header(), None, r.seq(), r.qual())
        };
        result.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}

fn main() {
    let matches = Command::new("subsmp")
        .about("Draw a random subset of records, or of reads by abundance.")
        .long_about(
            "Draw a random subset of records, or of reads by abundance.\n\n\
             With a seed, the draw is the one vsearch --fastx_subsample --randseed \
             gives; without one, or with 0, the seed comes from the system. The records \
             drawn are written in input order, unchanged, or with --by-size with the \
             number of their reads drawn as abundance (;size=). The whole input is read \
             and checked before anything is written.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA or FASTQ file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file of the records drawn, default to stdout"),
        )
        .arg(
            Arg::new("discarded")
                .long("discarded")
                .value_name("FILE")
                .help("Output file of the records, or reads, not drawn"),
        )
        .arg(
            Arg::new("count")
                .long("count")
                .value_name("N")
                .value_parser(value_parser!(u64).range(1..))
                .help("Draw this many records, or reads with --by-size"),
        )
        .arg(
            Arg::new("pct")
                .long("pct")
                .value_name("PCT")
                .value_parser(percentage)
                .help("Draw this percentage of records, or of reads with --by-size, rounded down"),
        )
        .group(ArgGroup::new("size").args(["count", "pct"]).required(true))
        .arg(
            Arg::new("by-size")
                .long("by-size")
                .action(ArgAction::SetTrue)
                .help(
                    "Count each record as its abundance (;size=), and write the reads drawn as it",
                ),
        )
        .arg(
            Arg::new("allow-fewer")
                .long("allow-fewer")
                .action(ArgAction::SetTrue)
                .help("Draw everything when asked for more than there is, instead of failing"),
        )
        .arg(
            Arg::new("seed")
                .long("seed")
                .value_name("SEED")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i64))
                .help("Seed of the random numbers, default to one from the system"),
        )
        .get_matches();
    let by_size = matches.get_flag("by-size");

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));
    let is_fasta = format == Some(Format::Fasta);

    // every record is read, and its abundance checked, before any is written
    let mut records = Vec::new();
    loop {
        let mut record = Record::new();
        if let Err(e) = reader.read(&mut record) {
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        let size = checked_abundance(record.header())
            .unwrap_or_else(|e| fail(&format!("{e} in read {}", record.header())));
        records.push((record, size));
    }

    let masses: Vec<u64> = records
        .iter()
        .map(|&(_, size)| if by_size { size } else { 1 })
        .collect();
    let total: u64 = masses.iter().sum();
    let wanted = match matches.get_one::<u64>("count") {
        Some(&count) => count,
        None => (total as f64 * matches.get_one::<f64>("pct").unwrap() / 100.0).floor() as u64,
    };
    if wanted > total && !matches.get_flag("allow-fewer") {
        fail(&format!(
            "Cannot subsample {wanted} reads from a sample of {total} (use --allow-fewer to keep them all)"
        ));
    }
    let mut generator = Mt19937_64::new(seed(matches.get_one::<i64>("seed").copied()));
    let drawn = draw(&masses, wanted.min(total), &mut generator);

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(io::stdout()),
    };
    write_counts(
        &mut Writer::new(fout, is_fasta, 0),
        &records,
        &drawn,
        by_size,
    );
    if let Some(f) = matches.get_one::<String>("discarded") {
        let rest: Vec<u64> = masses.iter().zip(&drawn).map(|(m, d)| m - d).collect();
        write_counts(
            &mut Writer::new(create(f), is_fasta, 0),
            &records,
            &rest,
            by_size,
        );
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
