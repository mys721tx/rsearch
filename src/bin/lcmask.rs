use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, ArgAction, Command};

use rsearch::dust::{Style, dust, hard_mask_lowercase};
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

/// What is done to each sequence.
#[derive(Clone, Copy)]
enum Masking {
    Dust(Style),
    LowercaseToN,
    Clear,
}

impl Masking {
    fn apply(self, seq: &mut [u8]) {
        match self {
            Masking::Dust(style) => dust(seq, style),
            Masking::LowercaseToN => hard_mask_lowercase(seq),
            Masking::Clear => seq.make_ascii_uppercase(),
        }
    }
}

fn main() {
    let matches = Command::new("lcmask")
        .about("Mask low-complexity regions with DUST.")
        .long_about(
            "Mask low-complexity regions with DUST.\n\n\
             By default the sequence is upper-cased and each low-complexity region \
             lower-cased (soft masking). With --hard each region becomes N instead and \
             the rest keeps its case; with --hard --no-dust only the lowercase residues \
             already there become N. --clear upper-cases every residue, unmasking it. \
             Headers and FASTQ scores are written unchanged, in the input format, FASTA \
             unwrapped.",
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
            Arg::new("hard")
                .long("hard")
                .action(ArgAction::SetTrue)
                .help("Mask with N instead of lowercase"),
        )
        .arg(
            Arg::new("no-dust")
                .long("no-dust")
                .action(ArgAction::SetTrue)
                .requires("hard")
                .help("Skip DUST: hard-mask only the lowercase residues already there"),
        )
        .arg(
            Arg::new("clear")
                .long("clear")
                .action(ArgAction::SetTrue)
                .conflicts_with_all(["hard", "no-dust"])
                .help("Upper-case every residue, clearing any masking"),
        )
        .get_matches();

    let masking = if matches.get_flag("clear") {
        Masking::Clear
    } else if matches.get_flag("no-dust") {
        Masking::LowercaseToN
    } else if matches.get_flag("hard") {
        Masking::Dust(Style::Hard)
    } else {
        Masking::Dust(Style::Soft)
    };

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
    let fail_after = |writer: &mut Writer<Box<dyn Write>>, message: &str| -> ! {
        // the reads before an error are written
        writer.flush().unwrap();
        fail(message);
    };

    let mut record = Record::new();
    let mut seq = Vec::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail_after(&mut writer, &e.to_string());
        }
        if record.is_empty() {
            break;
        }
        let header = record.header();
        // vsearch reads, and checks, every abundance
        if let Err(e) = checked_abundance(header) {
            fail_after(&mut writer, &format!("{e} in read {header}"));
        }

        seq.clear();
        seq.extend_from_slice(record.seq());
        masking.apply(&mut seq);
        writer
            .write(header, None, &seq, record.qual())
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }

    writer
        .flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
