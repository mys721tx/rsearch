use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
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
    let matches = Command::new("topn")
        .about("Keep the first N records.")
        .long_about(
            "Keep the first N records.\n\n\
             The whole input is read and checked, as vsearch does, so invalid input \
             after the first N records still fails. Records are written unchanged, in \
             the input format, FASTA unwrapped.",
        )
        .arg(
            Arg::new("count")
                .value_name("N")
                .required(true)
                .value_parser(value_parser!(u64).range(1..))
                .help("Number of records to keep"),
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
        .get_matches();
    let count = *matches.get_one::<u64>("count").unwrap();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, format == Some(Format::Fasta), 0);

    let mut record = Record::new();
    let mut written = 0;
    loop {
        if let Err(e) = reader.read(&mut record) {
            // the records before an error are written, as vsearch does
            writer.flush().unwrap();
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        // every record is read and checked; only the first are written
        if written < count {
            writer
                .write(record.header(), None, record.seq(), record.qual())
                .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
            written += 1;
        }
    }

    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
