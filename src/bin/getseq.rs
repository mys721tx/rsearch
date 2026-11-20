use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, Command};

use rsearch::hits::UcRecord;

use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::checked_abundance;
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

/// Selection by a UC file of hits: of the queries with hits, matched to
/// the records in order by label, or of the database records hit, by
/// their number.
enum Hits {
    Queries {
        lines: std::io::Lines<BufReader<File>>,
        pending: Option<String>,
    },
    Targets {
        hit: std::collections::HashSet<u64>,
        last: Option<u64>,
        seen: u64,
    },
}

impl Hits {
    fn new(reader: BufReader<File>, targets: bool) -> Hits {
        if !targets {
            return Hits::Queries {
                lines: reader.lines(),
                pending: None,
            };
        }
        let mut hit = std::collections::HashSet::new();
        for line in reader.lines() {
            let line = line.unwrap_or_else(|e| fail(&format!("cannot read the hits: {e}")));
            let record = UcRecord::parse(&line).unwrap_or_else(|e| fail(&e));
            if let Some(h) = record.hit {
                hit.insert(h.target);
            }
        }
        let last = hit.iter().max().copied();
        Hits::Targets { hit, last, seen: 0 }
    }

    /// whether the next record, of this header, is selected.
    fn select(&mut self, header: &str) -> Result<bool, String> {
        match self {
            Hits::Targets { hit, seen, .. } => {
                *seen += 1;
                Ok(hit.contains(&(*seen - 1)))
            }
            Hits::Queries { lines, pending } => {
                let prefix = format!("{header}\t");
                let mut selected = false;
                loop {
                    let line = match pending.take() {
                        Some(line) => line,
                        None => match lines.next() {
                            None => break,
                            Some(line) => line.map_err(|e| format!("cannot read the hits: {e}"))?,
                        },
                    };
                    let record = UcRecord::parse(&line)?;
                    if !record.labels.starts_with(&prefix) {
                        *pending = Some(line);
                        break;
                    }
                    if record.hit.is_some() {
                        selected = true;
                    } else {
                        // an N record is a query of its own
                        break;
                    }
                }
                Ok(selected)
            }
        }
    }

    /// fails if hits are left for queries or targets the input lacks.
    fn finish(self) -> Result<(), String> {
        match self {
            Hits::Targets { last, seen, .. } => match last {
                Some(t) if t >= seen => Err(format!("hit to target {t}, beyond the input")),
                _ => Ok(()),
            },
            Hits::Queries { mut lines, pending } => {
                match pending.or_else(|| lines.next().and_then(Result::ok)) {
                    Some(line) => Err(format!("hits of a query not in the input: {line}")),
                    None => Ok(()),
                }
            }
        }
    }
}

fn main() {
    let text = |name: &'static str, value: &'static str, help: &'static str| {
        Arg::new(name).long(name).value_name(value).help(help)
    };
    let mut command = Command::new("getseq")
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
        .arg(text(
            "hits",
            "FILE",
            "Select the queries with hits in this UC file, matched in order",
        ))
        .arg(
            Arg::new("targets")
                .long("targets")
                .action(ArgAction::SetTrue)
                .help("With --hits, select the database records hit, by their number"),
        )
        .group(
            ArgGroup::new("selection")
                .args(["label", "labels", "word", "words", "hits"])
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
        ));
    let matches = command.get_matches_mut();
    if matches.get_flag("targets") && matches.get_one::<String>("hits").is_none() {
        command
            .error(
                clap::error::ErrorKind::MissingRequiredArgument,
                "--targets needs --hits",
            )
            .exit();
    }

    let string = |name: &str| matches.get_one::<String>(name).map(String::as_str);
    let substr = matches.get_flag("substr");
    let mut hits = string("hits").map(|path| {
        let file = File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}")));
        Hits::new(BufReader::new(file), matches.get_flag("targets"))
    });
    let selector = if hits.is_some() {
        Selector::labels(std::iter::empty::<&str>(), false)
    } else if let Some(label) = string("label") {
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
        // vsearch reads, and checks, the abundance of every record
        if let Err(e) = checked_abundance(record.header()) {
            flush(&mut writer, &mut discarded);
            fail(&format!("{e} in read {}", record.header()));
        }
        let selected = match &mut hits {
            Some(hits) => hits.select(record.header()).unwrap_or_else(|e| {
                flush(&mut writer, &mut discarded);
                fail(&e)
            }),
            None => selector.matches(record.header()),
        };
        let out = if selected != invert {
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
    if let Some(hits) = hits {
        hits.finish().unwrap_or_else(|e| fail(&e));
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
