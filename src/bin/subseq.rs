use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};
use std::ops::Range;

use clap::{Arg, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::checked_abundance;
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

/// the bases from `start` to `end` of a read of `length` bases, both
/// counted from 1 and included; empty when the read ends before `start`.
fn region(start: u64, end: u64, length: usize) -> Range<usize> {
    let end = usize::try_from(end).unwrap_or(usize::MAX).min(length);
    let start = usize::try_from(start).unwrap_or(usize::MAX);
    if end < start { 0..0 } else { start - 1..end }
}

fn main() {
    let position = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("POS")
            .value_parser(value_parser!(u64).range(1..))
            .help(help)
    };
    let mut command = Command::new("subseq")
        .about("Keep one region of every read.")
        .long_about(
            "Keep one region of every read.\n\n\
             The region runs from --start to --end, both counted from 1 and included; \
             a read that ends before --start is written empty. Headers are written \
             unchanged, in the input format, FASTA unwrapped.",
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
        .arg(position("start", "First base kept, default to 1"))
        .arg(position(
            "end",
            "Last base kept, default to the end of the read",
        ));
    let matches = command.get_matches_mut();
    let start = matches.get_one::<u64>("start").copied().unwrap_or(1);
    let end = matches.get_one::<u64>("end").copied().unwrap_or(u64::MAX);
    if start > end {
        command
            .error(
                clap::error::ErrorKind::ValueValidation,
                "--start cannot be greater than --end",
            )
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
            break;
        }
        let header = record.header();
        // vsearch reads, and checks, the abundance of every record it writes
        if let Err(e) = checked_abundance(header) {
            fail_after(&mut writer, &format!("{e} in read {header}"));
        }
        let kept = region(start, end, record.seq().len());
        let qual = if is_fasta {
            &[][..]
        } else {
            &record.qual()[kept.clone()]
        };
        writer
            .write(header, None, &record.seq()[kept], qual)
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }

    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
