use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command};

use rsearch::cluster::reverse_complement_keeping_case;
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

fn open(path: &str) -> Box<dyn Read> {
    Box::new(File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}"))))
}

/// parses a pad, which must be ASCII.
fn ascii(text: &str) -> Result<String, String> {
    if text.is_ascii() {
        Ok(text.to_string())
    } else {
        Err("expected ASCII characters only".to_string())
    }
}

fn main() {
    let mut command = Command::new("joinpa")
        .about("Join read pairs with a pad.")
        .long_about(
            "Join read pairs with a pad.\n\n\
             Each forward read is followed by the pad and the reverse complement of the \
             read at the same position in --reverse, keeping its case; the scores are \
             those of the forward read, the pad and the reversed mate. The forward \
             header is kept whole. Both files are FASTQ, with as many reads each.",
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
                .help("Output FASTQ file of joined reads, default to stdout"),
        )
        .arg(
            Arg::new("reverse")
                .long("reverse")
                .value_name("FILE")
                .required(true)
                .help("Input FASTQ file of reverse reads"),
        )
        .arg(
            Arg::new("pad")
                .long("pad")
                .value_name("SEQ")
                .value_parser(ascii)
                .default_value("NNNNNNNN")
                .help("Sequence put between the mates"),
        )
        .arg(
            Arg::new("pad-qual")
                .long("pad-qual")
                .value_name("QUAL")
                .value_parser(ascii)
                .help("Scores of the pad, default to Q40 for each pad symbol"),
        )
        .arg(
            Arg::new("ascii")
                .long("ascii")
                .value_name("OFFSET")
                .value_parser(["33", "64"])
                .default_value("33")
                .help("ASCII offset of the default pad scores"),
        );
    let matches = command.get_matches_mut();

    let pad = matches
        .get_one::<String>("pad")
        .map_or("NNNNNNNN", String::as_str)
        .as_bytes()
        .to_vec();
    let pad_qual = match matches.get_one::<String>("pad-qual") {
        Some(qual) => qual.as_bytes().to_vec(),
        // Q40 at the offset
        None => {
            let symbol = match matches.get_one::<String>("ascii").map(String::as_str) {
                Some("64") => b'h',
                _ => b'I',
            };
            vec![symbol; pad.len()]
        }
    };
    if pad.len() != pad_qual.len() {
        command
            .error(
                clap::error::ErrorKind::ValueValidation,
                "--pad and --pad-qual must be equally long",
            )
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

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, false, 0);
    let fail_after = |writer: &mut Writer<Box<dyn Write>>, message: &str| -> ! {
        // the pairs before an error are written, as vsearch does
        let _ = writer.flush();
        fail(message);
    };

    let (mut record, mut mate) = (Record::new(), Record::new());
    let (mut seq, mut qual, mut complement) = (Vec::new(), Vec::new(), Vec::new());
    loop {
        if let Err(e) = forward.read(&mut record) {
            fail_after(&mut writer, &e.to_string());
        }
        if let Err(e) = reverse.read(&mut mate) {
            fail_after(&mut writer, &e.to_string());
        }
        match (record.is_empty(), mate.is_empty()) {
            (true, true) => break,
            (false, true) => fail_after(&mut writer, "More forward reads than reverse reads"),
            (true, false) => fail_after(&mut writer, "More reverse reads than forward reads"),
            (false, false) => {}
        }
        let header = record.header();
        // vsearch reads, and checks, the abundance of forward reads
        if let Err(e) = checked_abundance(header) {
            fail_after(&mut writer, &format!("{e} in read {header}"));
        }

        reverse_complement_keeping_case(&mut complement, mate.seq());
        seq.clear();
        seq.extend_from_slice(record.seq());
        seq.extend_from_slice(&pad);
        seq.extend_from_slice(&complement);
        qual.clear();
        qual.extend_from_slice(record.qual());
        qual.extend_from_slice(&pad_qual);
        qual.extend(mate.qual().iter().rev());
        writer
            .write(header, None, &seq, &qual)
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }

    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    for warning in [forward.stripped_warning(), reverse.stripped_warning()]
        .into_iter()
        .flatten()
    {
        warn(&warning);
    }
}
