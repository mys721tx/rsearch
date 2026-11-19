use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::accept::SizeFilters;
use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Reader, Record};
use rsearch::hits::{write_hit, write_no_hit};
use rsearch::label::checked_abundance;

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

/// a sequence as vsearch compares it: upper case, U as T.
fn normalize(out: &mut Vec<u8>, seq: &[u8]) {
    out.clear();
    out.extend(seq.iter().map(|&c| match c.to_ascii_uppercase() {
        b'U' => b'T',
        c => c,
    }));
}

/// The database: labels, abundances, and the numbers of the sequences of
/// each normalized sequence, in order.
struct Database {
    labels: Vec<String>,
    sizes: Vec<i64>,
    by_sequence: HashMap<Vec<u8>, Vec<u32>>,
}

fn database(path: &str) -> Database {
    let file = File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}")));
    let mut reader = Reader::new(BufReader::new(file));
    let mut db = Database {
        labels: Vec::new(),
        sizes: Vec::new(),
        by_sequence: HashMap::new(),
    };
    let mut record = Record::new();
    let mut normalized = Vec::new();
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            break;
        }
        let size = checked_abundance(record.header())
            .unwrap_or_else(|e| fail(&format!("{e} in read {} of the database", record.header())));
        normalize(&mut normalized, record.seq());
        db.by_sequence
            .entry(normalized.clone())
            .or_default()
            .push(db.labels.len() as u32);
        db.labels.push(record.header().to_string());
        db.sizes.push(size as i64);
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    db
}

fn main() {
    let ratio = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("RATIO")
            .value_parser(value_parser!(f64))
            .help(help)
    };
    let size = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("SIZE")
            .value_parser(value_parser!(i64).range(0..))
            .help(help)
    };
    let matches = Command::new("exsrch")
        .about("Find exact full-length matches in a database by hashing; write hits.")
        .long_about(
            "Find exact full-length matches in a database by hashing; write hits.\n\n\
             As vsearch --search_exact: the hits of a query are the database sequences \
             identical to it, ignoring case and taking U for T, by target number. They \
             are written in vsearch's UC format, every hit as with --uc_allhits, and an \
             N record for a query without hits.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA or FASTQ file of queries, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output UC file of hits, default to stdout"),
        )
        .arg(
            Arg::new("db")
                .long("db")
                .value_name("FILE")
                .required(true)
                .help("Database: FASTA or FASTQ file"),
        )
        .arg(
            Arg::new("strand")
                .long("strand")
                .value_name("STRAND")
                .value_parser(["plus", "both"])
                .default_value("plus")
                .help("Search the plus strand, or both"),
        )
        .arg(
            Arg::new("maxhits")
                .long("maxhits")
                .value_name("N")
                .value_parser(value_parser!(u64).range(1..))
                .help("Write at most this many hits of a query, default to all"),
        )
        .arg(
            Arg::new("self")
                .long("self")
                .action(ArgAction::SetTrue)
                .help("Skip targets with the query's label"),
        )
        .arg(size("max-qsize", "Skip queries of an abundance above this"))
        .arg(size("min-tsize", "Skip targets of an abundance below this"))
        .arg(ratio(
            "min-size-ratio",
            "Skip targets unless the query's abundance is at least this times theirs",
        ))
        .arg(ratio(
            "max-size-ratio",
            "Skip targets unless the query's abundance is at most this times theirs",
        ))
        .get_matches();
    let both_strands = matches.get_one::<String>("strand").unwrap() == "both";
    let max_hits = matches
        .get_one::<u64>("maxhits")
        .copied()
        .unwrap_or(u64::MAX);
    let skip_self = matches.get_flag("self");
    let mut filters = SizeFilters::default();
    if let Some(&v) = matches.get_one::<i64>("max-qsize") {
        filters.max_query_size = v;
    }
    if let Some(&v) = matches.get_one::<i64>("min-tsize") {
        filters.min_target_size = v;
    }
    if let Some(&v) = matches.get_one::<f64>("min-size-ratio") {
        filters.min_size_ratio = v;
    }
    if let Some(&v) = matches.get_one::<f64>("max-size-ratio") {
        filters.max_size_ratio = v;
    }

    // the database is read whole, then the output created and the queries
    // read, as vsearch does
    let db = database(matches.get_one::<String>("db").unwrap());
    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    let (mut normalized, mut rc) = (Vec::new(), Vec::new());
    let mut hits: Vec<(u32, bool)> = Vec::new();
    loop {
        let size = match reader.read(&mut record) {
            Err(e) => Err(e.to_string()),
            Ok(()) if record.is_empty() => break,
            Ok(()) => checked_abundance(record.header())
                .map_err(|e| format!("{e} in read {}", record.header())),
        };
        let size = size.unwrap_or_else(|message| {
            // the hits of the queries before an error are written
            let _ = out.flush();
            fail(&message)
        }) as i64;
        hits.clear();
        let header = record.header();
        let mut search = |seq: &[u8], minus: bool, hits: &mut Vec<(u32, bool)>| {
            normalize(&mut normalized, seq);
            for &target in db.by_sequence.get(&normalized).into_iter().flatten() {
                let t = target as usize;
                if filters.accepts(size, db.sizes[t]) && !(skip_self && db.labels[t] == header) {
                    hits.push((target, minus));
                }
            }
        };
        search(record.seq(), false, &mut hits);
        if both_strands {
            reverse_complement_keeping_case(&mut rc, record.seq());
            search(&rc, true, &mut hits);
        }
        // by target, the plus strand first
        hits.sort_by_key(|&(target, minus)| (target, minus));
        let written = if hits.is_empty() {
            write_no_hit(&mut out, header.as_bytes())
        } else {
            hits.iter()
                .take(usize::try_from(max_hits).unwrap_or(usize::MAX))
                .try_for_each(|&(target, minus)| {
                    write_hit(
                        &mut out,
                        u64::from(target),
                        record.seq().len(),
                        100.0,
                        minus,
                        b"=",
                        header.as_bytes(),
                        db.labels[target as usize].as_bytes(),
                    )
                })
        };
        written.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    out.flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
}
