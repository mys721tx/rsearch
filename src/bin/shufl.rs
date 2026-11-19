use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::random::{Mt19937_64, seed, shuffle};
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

fn main() {
    let matches = Command::new("shufl")
        .about("Shuffle the order of records.")
        .long_about(
            "Shuffle the order of records.\n\n\
             With a seed, the order is the one vsearch --shuffle --randseed gives; without \
             one, or with 0, the seed comes from the system. The whole input is read and \
             checked before anything is written. Records are written unchanged, in the \
             input format, FASTA unwrapped.",
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
                .help("Output file in the input format, default to stdout"),
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

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));

    // every record is read, and checked, before any is written
    let mut records = Vec::new();
    loop {
        let mut record = Record::new();
        if let Err(e) = reader.read(&mut record) {
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        if let Err(e) = checked_abundance(record.header()) {
            fail(&format!("{e} in read {}", record.header()));
        }
        records.push(record);
    }

    let mut generator = Mt19937_64::new(seed(matches.get_one::<i64>("seed").copied()));
    shuffle(&mut records, &mut generator);

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, format == Some(Format::Fasta), 0);
    for r in &records {
        writer
            .write(r.header(), None, r.seq(), r.qual())
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
