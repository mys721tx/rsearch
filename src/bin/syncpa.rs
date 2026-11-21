use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, Command};

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

fn open(path: &str) -> Box<dyn Read> {
    Box::new(File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}"))))
}

fn create(path: &str) -> Box<dyn Write> {
    Box::new(File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))))
}

/// The key mates share, as vsearch derives it: the header up to its first
/// space or tab, less a final 1 or 2 after one of `separators`.
///
/// Casava 1.8 headers differ only after the space; older ones end in /1
/// and /2.
fn key<'a>(header: &'a [u8], separators: &[u8]) -> &'a [u8] {
    let end = header
        .iter()
        .position(|&c| c == b' ' || c == b'\t')
        .unwrap_or(header.len());
    let key = &header[..end];
    match key {
        [.., separator, b'1' | b'2'] if separators.contains(separator) => &key[..key.len() - 2],
        _ => key,
    }
}

/// The four outputs; all but the synced forward reads optional.
struct Outputs {
    forward: Writer<Box<dyn Write>>,
    reverse: Option<Writer<Box<dyn Write>>>,
    orphans: Option<Writer<Box<dyn Write>>>,
    rev_orphans: Option<Writer<Box<dyn Write>>>,
}

impl Outputs {
    fn all(&mut self) -> impl Iterator<Item = &mut Writer<Box<dyn Write>>> {
        std::iter::once(&mut self.forward).chain(
            [&mut self.reverse, &mut self.orphans, &mut self.rev_orphans]
                .into_iter()
                .flatten(),
        )
    }

    /// writes the reads so far, then fails, as vsearch leaves them.
    fn fail(&mut self, message: &str) -> ! {
        for writer in self.all() {
            let _ = writer.flush();
        }
        fail(message);
    }
}

/// writes a record to an output, if there is one.
fn write(writer: Option<&mut Writer<Box<dyn Write>>>, record: &Record) {
    if let Some(writer) = writer {
        writer
            .write(record.header(), None, record.seq(), record.qual())
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
}

fn main() {
    let file = |name: &'static str, help: &'static str| {
        Arg::new(name).long(name).value_name("FILE").help(help)
    };
    let matches = Command::new("syncpa")
        .about("Pair the reads of two files by label.")
        .long_about(
            "Pair the reads of two files by label.\n\n\
             Mates share their header up to its first space or tab, less a final 1 or 2 \
             after one of --separators. Forward reads with a mate are written in their \
             order, and their mates in the same order to --rev-out; reads without a mate \
             go to --orphans and --rev-orphans. Reads are written as read, in the input \
             format, FASTA unwrapped.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input file of forward reads, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file of forward reads with a mate, default to stdout"),
        )
        .arg(file("reverse", "Input file of reverse reads").required(true))
        .arg(file(
            "rev-out",
            "Output file of the mates, in the order of the forward reads",
        ))
        .arg(file(
            "orphans",
            "Output file of forward reads without a mate",
        ))
        .arg(file(
            "rev-orphans",
            "Output file of reverse reads without a mate",
        ))
        .arg(
            Arg::new("separators")
                .long("separators")
                .value_name("CHARS")
                .default_value("/")
                .help("Characters that may come before a final 1 or 2"),
        )
        .get_matches();
    let separators = matches
        .get_one::<String>("separators")
        .map_or(&b"/"[..], |s| s.as_bytes());

    let fin = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => open(f),
        None => Box::new(io::stdin()),
    };
    let mut forward = Reader::new(BufReader::new(fin));
    let format = forward.format().unwrap_or_else(|e| fail(&e.to_string()));
    let mut reverse = Reader::new(BufReader::new(open(
        matches.get_one::<String>("reverse").unwrap(),
    )));
    let rev_format = reverse.format().unwrap_or_else(|e| fail(&e.to_string()));
    if format.is_some() && rev_format.is_some() && format != rev_format {
        fail("Forward and reverse files must both be FASTA or both FASTQ");
    }
    // the format of the file with reads
    let is_fasta = format.or(rev_format) == Some(Format::Fasta);

    // the reverse reads, in order, and the position of each key
    let mut mates = Vec::new();
    let mut positions: HashMap<Vec<u8>, usize> = HashMap::new();
    loop {
        let mut record = Record::new();
        if let Err(e) = reverse.read(&mut record) {
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        let key = key(record.header().as_bytes(), separators).to_vec();
        if positions.insert(key, mates.len()).is_some() {
            fail("Duplicate read label in reverse file");
        }
        mates.push(record);
    }

    let writer = |out: Box<dyn Write>| Writer::new(out, is_fasta, 0);
    let output = |name: &str| matches.get_one::<String>(name).map(|f| writer(create(f)));
    let mut outputs = Outputs {
        forward: writer(match matches.get_one::<String>("out").map(String::as_str) {
            Some(f) => create(f),
            None => Box::new(io::stdout()),
        }),
        reverse: output("rev-out"),
        orphans: output("orphans"),
        rev_orphans: output("rev-orphans"),
    };

    let mut paired = vec![false; mates.len()];
    let mut record = Record::new();
    loop {
        if let Err(e) = forward.read(&mut record) {
            outputs.fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        match positions.get(key(record.header().as_bytes(), separators)) {
            None => write(outputs.orphans.as_mut(), &record),
            Some(&position) => {
                // two forward reads for one mate cannot be paired
                if paired[position] {
                    outputs.fail("Duplicate read label in forward file");
                }
                paired[position] = true;
                write(Some(&mut outputs.forward), &record);
                write(outputs.reverse.as_mut(), &mates[position]);
            }
        }
    }
    for (mate, _) in mates.iter().zip(&paired).filter(|(_, paired)| !**paired) {
        write(outputs.rev_orphans.as_mut(), mate);
    }

    for writer in outputs.all() {
        writer
            .flush()
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    for reader_warning in [forward.stripped_warning(), reverse.stripped_warning()]
        .into_iter()
        .flatten()
    {
        warn(&reader_warning);
    }
}
