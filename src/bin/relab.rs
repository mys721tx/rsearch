use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, ArgAction, Command};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::filter::Filter;
use rsearch::label::{Label, checked_abundance};
use rsearch::options;
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
    let command = Command::new("relab")
        .about("Rewrite sequence labels: relabel, add or remove annotations.")
        .long_about(
            "Rewrite sequence labels: relabel, add or remove annotations.\n\n\
             Labels are replaced by a prefix and the record's number, the sequence, \
             or the SHA1 or MD5 digest of the normalized sequence; or kept, with \
             ;size=, ;ee= and ;length= annotations removed. Then the suffix, \
             ;sample=, ;size=, ;ee= and ;length= are appended, in that order. \
             Records are written in the input format, FASTA unwrapped.",
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
                .help("Output file, default to stdout"),
        )
        .arg(
            Arg::new("sizeout")
                .long("sizeout")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with abundance (;size=)"),
        )
        .arg(
            Arg::new("ee")
                .long("ee")
                .action(ArgAction::SetTrue)
                .help("Annotate output labels with expected errors (;ee=), FASTQ only"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .requires("ee")
                .default_value("33")
                .help("ASCII offset of quality scores for --ee"),
        );
    let mut command = options::label_args(command);
    let matches = command.get_matches_mut();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));
    let is_fasta = format == Some(Format::Fasta);

    let label = Label {
        sizeout: matches.get_flag("sizeout"),
        eeout: matches.get_flag("ee"),
        ..options::label_from(&matches)
    };
    if is_fasta && label.eeout {
        command
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                "--ee cannot be used with FASTA input",
            )
            .exit();
    }
    let offset: u8 = matches
        .get_one::<String>("ascii")
        .map_or(33, |a| a.parse().unwrap());
    // expected errors of the whole read, scores from 0 to the highest
    // printable
    let scores = Filter {
        min_len: 0,
        offset,
        max_score: 126 - i16::from(offset),
        ..Filter::default()
    };

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, is_fasta, 0);
    let fail_after = |writer: &mut Writer<Box<dyn Write>>, message: &str| -> ! {
        // the records before an error are written, as vsearch does
        writer.flush().unwrap();
        fail(message);
    };

    let mut record = Record::new();
    let mut ordinal = 0;
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
        ordinal += 1;

        let header = record.header();
        // vsearch reads, and checks, every abundance
        let size = checked_abundance(header)
            .unwrap_or_else(|e| fail_after(&mut writer, &format!("{e} in read {header}")));
        let ee = if label.eeout {
            match scores.apply(record.seq(), Some(record.qual()), size) {
                Ok(outcome) => outcome.ee,
                Err(e) => fail_after(&mut writer, &format!("{e} in read {}", record.id())),
            }
        } else {
            None
        };

        let result = if label.is_identity() {
            writer.write(header, None, record.seq(), record.qual())
        } else {
            let header = label.format(header, record.seq(), size, ordinal, ee);
            writer.write(&header, None, record.seq(), record.qual())
        };
        result.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }

    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
