use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, Command};

use rsearch::fastx::{Format, Reader, Record};
use rsearch::select::{Selector, labels_of};
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

fn create(path: &str) -> Box<dyn Write> {
    Box::new(File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))))
}

/// reads a file of labels or words, one per line.
fn read_list(path: &str) -> Vec<Vec<u8>> {
    let text = std::fs::read(path)
        .unwrap_or_else(|e| fail(&format!("Unable to open labels file ({path}): {e}")));
    labels_of(&text).into_iter().map(<[u8]>::to_vec).collect()
}

fn main() {
    let text = |name: &'static str, value: &'static str, help: &'static str| {
        Arg::new(name).long(name).value_name(value).help(help)
    };
    let matches = Command::new("getseq")
        .about("Select records by their headers.")
        .long_about(
            "Select records by their headers.\n\n\
             A label matches the whole header, ignoring case, or with --substr any part \
             of it. A word matches between characters that are not letters or digits, \
             case included, or with --field the value of NAME=word between semicolons. \
             Records are written unchanged, in the input format, FASTA unwrapped.",
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
                .help("Output file of selected records, default to stdout"),
        )
        .arg(text("label", "LABEL", "Select records with this header"))
        .arg(text(
            "labels",
            "FILE",
            "Select records with a header listed in FILE",
        ))
        .arg(text(
            "word",
            "WORD",
            "Select records whose header holds this word",
        ))
        .arg(text(
            "words",
            "FILE",
            "Select records whose header holds a word listed in FILE",
        ))
        .group(
            ArgGroup::new("selection")
                .args(["label", "labels", "word", "words"])
                .required(true),
        )
        .arg(
            Arg::new("substr")
                .long("substr")
                .action(ArgAction::SetTrue)
                .conflicts_with_all(["word", "words"])
                .help("Match labels as any part of the header"),
        )
        .arg(
            text(
                "field",
                "NAME",
                "Match words as the value of NAME= between semicolons",
            )
            .conflicts_with_all(["label", "labels"]),
        )
        .arg(
            Arg::new("invert")
                .long("invert")
                .action(ArgAction::SetTrue)
                .help("Write the records that are not selected instead"),
        )
        .arg(text(
            "discarded",
            "FILE",
            "Output file of the records not written to the output",
        ))
        .get_matches();

    let string = |name: &str| matches.get_one::<String>(name).map(String::as_str);
    let substr = matches.get_flag("substr");
    let selector = if let Some(label) = string("label") {
        Selector::labels([label], substr)
    } else if let Some(path) = string("labels") {
        Selector::labels(read_list(path), substr)
    } else if let Some(word) = string("word") {
        Selector::words([word], string("field"))
    } else {
        Selector::words(read_list(string("words").unwrap()), string("field"))
    };
    let invert = matches.get_flag("invert");

    let fin: Box<dyn Read> = match string("in") {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let format = reader.format().unwrap_or_else(|e| fail(&e.to_string()));
    let is_fasta = format == Some(Format::Fasta);

    let fout: Box<dyn Write> = match string("out") {
        Some(f) => create(f),
        None => Box::new(io::stdout()),
    };
    let mut writer = Writer::new(fout, is_fasta, 0);
    let mut discarded = string("discarded").map(|f| Writer::new(create(f), is_fasta, 0));
    let flush = |writer: &mut Writer<Box<dyn Write>>,
                 discarded: &mut Option<Writer<Box<dyn Write>>>| {
        for w in [Some(writer), discarded.as_mut()].into_iter().flatten() {
            w.flush()
                .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
        }
    };

    let mut record = Record::new();
    loop {
        if let Err(e) = reader.read(&mut record) {
            // the records before an error are written, as vsearch does
            flush(&mut writer, &mut discarded);
            fail(&e.to_string());
        }
        if record.is_empty() {
            break;
        }
        let out = if selector.matches(record.header()) != invert {
            Some(&mut writer)
        } else {
            discarded.as_mut()
        };
        if let Some(out) = out {
            out.write(record.header(), None, record.seq(), record.qual())
                .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
        }
    }

    flush(&mut writer, &mut discarded);
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
