use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::random::{SplitMix64, seed, substream};
use rsearch::scramble::Scrambler;
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
    let matches = Command::new("scrmbl")
        .about("Shuffle the bases within each read.")
        .long_about(
            "Shuffle the bases within each read.\n\n\
             With --kmer 1, the bases are shuffled; with a k from 2 to 9, they are \
             rearranged at random keeping the count of every word of up to k bases. \
             With a seed, each read is scrambled as vsearch --scramble --randseed \
             scrambles it; without one, or with 0, the seed comes from the system. \
             Scores and headers are written unchanged, in the input format, FASTA \
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
            Arg::new("kmer")
                .long("kmer")
                .value_name("K")
                .value_parser(value_parser!(u8).range(1..=9))
                .help("Keep the counts of words of up to K bases, 1 to 9, default to 1"),
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
    let k = usize::from(matches.get_one::<u8>("kmer").copied().unwrap_or(1));
    let base_seed = seed(matches.get_one::<i64>("seed").copied());

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
        // the reads before an error are written, as vsearch does
        writer.flush().unwrap();
        fail(message);
    };

    let mut scrambler = Scrambler::default();
    let mut record = Record::new();
    let mut seq = Vec::new();
    let mut count = 0;
    loop {
        if let Err(e) = reader.read(&mut record) {
            fail_after(&mut writer, &e.to_string());
        }
        if record.is_empty() {
            break;
        }
        count += 1;
        let header = record.header();
        if let Err(e) = checked_abundance(header) {
            fail_after(&mut writer, &format!("{e} in read {header}"));
        }
        // each read's numbers depend only on the seed and its number
        let mut generator = SplitMix64::new(substream(base_seed, count));
        seq.clear();
        seq.extend_from_slice(record.seq());
        scrambler.scramble(&mut seq, k, &mut generator);
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
