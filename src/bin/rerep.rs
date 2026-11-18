use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, ArgAction, Command};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::{Label, checked_abundance, size_annotation};
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
    let matches = Command::new("rerep")
        .about("Write each read as many times as its abundance (;size=).")
        .long_about(
            "Write each read as many times as its abundance (;size=).\n\n\
             Each copy is annotated with an abundance of 1, replacing the read's; a \
             read without an abundance is written once. Records are written in the \
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
            Arg::new("no-sizeout")
                .long("no-sizeout")
                .action(ArgAction::SetTrue)
                .help("Remove abundances (;size=) instead of writing ;size=1"),
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
    let is_fasta = format == Some(Format::Fasta);

    // the copies' label: the read's, with an abundance of 1 or none
    let label = if matches.get_flag("no-sizeout") {
        Label {
            xsize: true,
            ..Label::default()
        }
    } else {
        Label {
            sizeout: true,
            ..Label::default()
        }
    };

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
    let mut missing = false;
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail_after(&mut writer, &e.to_string());
        }
        if record.is_empty() {
            break;
        }

        let header = record.header();
        let copies = checked_abundance(header)
            .unwrap_or_else(|e| fail_after(&mut writer, &format!("{e} in read {header}")));
        // after a valid abundance, none means an absent one
        missing |= size_annotation(header).is_none();

        let copy = label.format(header, record.seq(), 1, 0, None);
        for _ in 0..copies {
            writer
                .write(&copy, None, record.seq(), record.qual())
                .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
        }
    }

    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if missing {
        warn(b"Missing abundance information for some input sequences, assumed 1");
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
