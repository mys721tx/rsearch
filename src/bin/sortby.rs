use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::cmp::Ordering;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, Command};

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

/// A read with its abundance.
struct Entry {
    record: Record,
    size: u64,
}

/// orders reads by decreasing abundance, then header, as vsearch's
/// `--sortbysize`. The sort is stable, and input order breaks the last
/// ties.
fn by_size(a: &Entry, b: &Entry) -> Ordering {
    b.size.cmp(&a.size).then_with(|| {
        a.record
            .header()
            .as_bytes()
            .cmp(b.record.header().as_bytes())
    })
}

/// orders reads by decreasing length, then as `by_size`, as vsearch's
/// `--sortbylength`.
fn by_length(a: &Entry, b: &Entry) -> Ordering {
    b.record
        .seq()
        .len()
        .cmp(&a.record.seq().len())
        .then_with(|| by_size(a, b))
}

fn main() {
    let matches = Command::new("sortby")
        .about("Sort reads by decreasing length or abundance.")
        .long_about(
            "Sort reads by decreasing length or abundance.\n\n\
             By length, ties are broken by decreasing abundance (;size=), then header, \
             then input order; by abundance, by header, then input order. The whole \
             input is read and checked before anything is written. Records are \
             written unchanged, in the input format, FASTA unwrapped.",
        )
        .arg(
            Arg::new("length")
                .long("length")
                .action(ArgAction::SetTrue)
                .help("Sort by decreasing length, as vsearch --sortbylength"),
        )
        .arg(
            Arg::new("size")
                .long("size")
                .action(ArgAction::SetTrue)
                .help("Sort by decreasing abundance (;size=), as vsearch --sortbysize"),
        )
        .group(ArgGroup::new("key").args(["length", "size"]).required(true))
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

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));

    // every read is read, and checked, before any is written
    let mut reads = Vec::new();
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
        reads.push(Entry { record, size });
    }

    // a stable sort, so that input order breaks the last ties
    if matches.get_flag("length") {
        reads.sort_by(by_length);
    } else {
        reads.sort_by(by_size);
    }

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, format == Some(Format::Fasta), 0);
    for Entry { record: r, .. } in &reads {
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
