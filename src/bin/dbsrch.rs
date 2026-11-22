use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Reader, Record};
use rsearch::hits::{write_hit, write_no_hit};
use rsearch::label::checked_abundance;
use rsearch::options::{hit_choice_args, search_args, search_options};
use rsearch::search::{Limits, MIN_WORD_MATCHES, Searcher, order, read_database};
use rsearch::threads::{Workers, threads, threads_args};
use rsearch::words::{MAX_WORD_LEN, MIN_WORD_LEN};

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
    let command = Command::new("dbsrch")
        .about("Search queries against a database; write hits.")
        .long_about(
            "Search queries against a database; write hits.\n\n\
             As vsearch --usearch_global: the database targets sharing the most words \
             with a query (lower case masked) are aligned globally with it, best first, \
             until --maxaccepts are accepted or --maxrejects rejected; the hits accepted \
             and the weak ones are written in UC format, by decreasing identity, then \
             target, and an N record for a query without hits.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA or FASTQ queries, default to stdin"),
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
                .value_name("DATABASE")
                .required(true)
                .help("Database: FASTA, FASTQ or UDB file"),
        )
        .arg(
            Arg::new("strand")
                .long("strand")
                .value_name("STRAND")
                .value_parser(["plus", "both"])
                .default_value("plus")
                .help("Search the plus strand of the queries, or both"),
        )
        .arg(
            Arg::new("maxaccepts")
                .long("maxaccepts")
                .value_name("N")
                .value_parser(value_parser!(u64))
                .default_value("1")
                .help("Stop after this many accepted hits a strand, 0 for no limit"),
        )
        .arg(
            Arg::new("maxrejects")
                .long("maxrejects")
                .value_name("N")
                .value_parser(value_parser!(u64))
                .default_value("32")
                .help("Stop after this many rejected targets a strand, 0 for no limit"),
        )
        .arg(
            Arg::new("word-len")
                .long("word-len")
                .value_name("LENGTH")
                .value_parser(
                    value_parser!(u32).range(i64::from(MIN_WORD_LEN)..=i64::from(MAX_WORD_LEN)),
                )
                .default_value("8")
                .help("Length of the words indexed (a UDB database's own wins)"),
        )
        .arg(
            Arg::new("min-word-matches")
                .long("min-word-matches")
                .value_name("N")
                .value_parser(value_parser!(u32))
                .help(
                    "Words a target shares with the query to be a candidate, default by --word-len",
                ),
        );
    let matches = threads_args(hit_choice_args(search_args(command)))
        .mut_arg("id", |a| a.required(true))
        .get_matches();
    let (scoring, criteria) = search_options(&matches);
    let iddef = criteria.iddef;
    let max_hits = matches
        .get_one::<u64>("maxhits")
        .map_or(usize::MAX, |&n| n as usize);
    let top_only = matches.get_flag("top-hits-only");
    let both = matches.get_one::<String>("strand").unwrap() == "both";
    let word_len = *matches.get_one::<u32>("word-len").unwrap();
    // vsearch takes the default from the word length asked for, before a
    // UDB file's own
    let min_word_matches = matches
        .get_one::<u32>("min-word-matches")
        .copied()
        .unwrap_or(MIN_WORD_MATCHES[word_len as usize]);
    let limits = Limits {
        max_accepts: *matches.get_one::<u64>("maxaccepts").unwrap() as usize,
        max_rejects: *matches.get_one::<u64>("maxrejects").unwrap() as usize,
        min_word_matches,
        top_hits: None,
    };

    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let db = read_database(
        matches.get_one::<String>("db").unwrap(),
        word_len,
        &mut |w| warn(w),
    )
    .unwrap_or_else(|e| fail(&e));
    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    // a query, or the error that ends the input. Nothing is read after an
    // error, and the hits of the queries before it are written first
    let mut failed = false;
    let next = || -> Option<Result<(String, i64, Vec<u8>), String>> {
        if failed {
            return None;
        }
        let query = match reader.read(&mut record) {
            Err(e) => Err(e.to_string()),
            Ok(()) if record.is_empty() => return None,
            Ok(()) => {
                let label = record.header();
                checked_abundance(label)
                    .map(|size| (label.to_string(), size as i64, record.seq().to_vec()))
                    .map_err(|e| format!("{e} in read {label}"))
            }
        };
        failed = query.is_err();
        Some(query)
    };
    // each worker's searcher and reverse strand
    let threads = threads(&matches);
    let workers = Workers::new(threads, || {
        (
            Searcher::new(&db, criteria.clone(), scoring, limits),
            Vec::new(),
        )
    });
    let search = |(searcher, minus): &mut (Searcher, Vec<u8>),
                  query: Result<(String, i64, Vec<u8>), String>|
     -> Result<Vec<u8>, String> {
        let (label, qsize, query) = query?;
        let mut hits = searcher.search(&db, &label, qsize, &query, false);
        if both {
            reverse_complement_keeping_case(minus, &query);
            hits.extend(searcher.search(&db, &label, qsize, minus, true));
        }
        order(&mut hits, iddef);
        hits.truncate(max_hits);
        if top_only && let Some(best) = hits.first().map(|h| h.stats.id[iddef]) {
            let n = hits
                .iter()
                .take_while(|h| h.stats.id[iddef] >= best)
                .count();
            hits.truncate(n);
        }
        let mut text = Vec::new();
        if hits.is_empty() {
            write_no_hit(&mut text, label.as_bytes()).unwrap();
        } else {
            for h in &hits {
                let perfect = h.stats.matches == h.stats.length;
                write_hit(
                    &mut text,
                    h.target as u64,
                    query.len(),
                    h.stats.id[iddef],
                    h.minus,
                    if perfect { b"=" } else { h.cigar.as_bytes() },
                    label.as_bytes(),
                    db.labels[h.target].as_bytes(),
                )
                .unwrap();
            }
        }
        Ok(text)
    };
    workers.stream(32 * threads, next, search, |text| match text {
        Ok(text) => out
            .write_all(&text)
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}"))),
        Err(message) => fail(&message),
    });
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    out.flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
