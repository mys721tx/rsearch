use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::fastx::{Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::udb;
use rsearch::words::{Index, MAX_WORD_LEN, MIN_WORD_LEN};

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
    let matches = Command::new("mkudb")
        .about("Build a UDB database: sequences and word index.")
        .long_about(
            "Build a UDB database: sequences and word index.\n\n\
             Writes vsearch's UDB file of the input: the sequences as given, headers \
             whole, and for each word of --word-len residues the sequences holding it. \
             A word holding lower case (masked) residues or symbols other than A, C, G, \
             T or U is not indexed; mask with lcmask first. FASTA or FASTQ in, the UDB \
             file out.",
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
                .help("Output UDB file, default to stdout"),
        )
        .arg(
            Arg::new("word-len")
                .long("word-len")
                .value_name("LEN")
                .default_value("8")
                .value_parser(
                    value_parser!(u32).range(i64::from(MIN_WORD_LEN)..=i64::from(MAX_WORD_LEN)),
                )
                .help("Length of the indexed words, 3 to 15"),
        )
        .get_matches();
    let word_len = *matches.get_one::<u32>("word-len").unwrap();

    // the output is created first, as vsearch does, so that a failure
    // leaves it empty
    let mut out: Box<dyn Write> = match matches.get_one::<String>("out") {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));

    // every record is read, and checked, before anything is written
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut record = Record::new();
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            break;
        }
        if let Err(e) = checked_abundance(record.header()) {
            fail(&format!("{e} in read {}", record.header()));
        }
        records.push((record.header().to_string(), record.seq().to_vec()));
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    if records.is_empty() {
        // vsearch reads no UDB file without sequences
        fail("no sequences: a UDB file needs at least one");
    }

    let index = Index::new(word_len, records.iter().map(|r| r.1.as_slice()));
    if index.without_words() > 0 {
        warn(
            format!(
                "{} of {} sequences give no word to index: they are shorter than the \
                 word length ({word_len}), or each of their words holds a lower-case \
                 residue or a symbol other than A, C, G, T or U",
                index.without_words(),
                records.len()
            )
            .as_bytes(),
        );
    }
    udb::write(&mut out, &index, &records)
        .and_then(|()| out.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
