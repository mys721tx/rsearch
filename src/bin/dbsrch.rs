use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Reader, Record};
use rsearch::hits::{write_hit, write_no_hit};
use rsearch::label::checked_abundance;
use rsearch::options::{search_args, search_options};
use rsearch::search::{Database, Limits, MIN_WORD_MATCHES, Searcher, order};
use rsearch::udb::{self, Header, Udb};
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

/// the abundance of a database record, which vsearch always reads.
fn size(header: &str) -> i64 {
    checked_abundance(header).unwrap_or_else(|e| fail(&format!("{e} in the database"))) as i64
}

/// reads the database: a UDB file, whose word length wins over
/// `word_len`, or a FASTA or FASTQ file indexed at `word_len`.
fn database(path: &str, word_len: u32) -> Database {
    let file = File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}")));
    let file_size = file
        .metadata()
        .ok()
        .filter(|m| m.is_file())
        .map(|m| m.len());
    let mut fin = BufReader::new(file);
    let start = fin
        .fill_buf()
        .unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
    if start.starts_with(udb::MAGIC) {
        let mut bytes = Vec::new();
        let file_size = match file_size {
            Some(s) => s,
            None => {
                fin.read_to_end(&mut bytes)
                    .unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
                bytes.len() as u64
            }
        };
        // vsearch warns as it reads the header
        let first = if bytes.is_empty() {
            fin.fill_buf().unwrap_or_default()
        } else {
            &bytes
        };
        if let Ok(header) = Header::parse(first, file_size)
            && header.word_len() != word_len
        {
            warn(
                format!(
                    "the word length of the UDB file, {}, is used",
                    header.word_len()
                )
                .as_bytes(),
            );
        }
        let mut list = Vec::new();
        let keep = |_: &[u32], _: u64, block: &[u32]| list.extend_from_slice(block);
        let read = if bytes.is_empty() {
            Udb::read(&mut fin, file_size, keep)
        } else {
            Udb::read(&mut &bytes[..], file_size, keep)
        };
        let udb = read.unwrap_or_else(|e| fail(&e.to_string()));
        let labels: Vec<String> = udb
            .header_starts
            .windows(2)
            .map(|pair| {
                String::from_utf8_lossy(&udb.headers[pair[0] as usize..pair[1] as usize - 1])
                    .into_owned()
            })
            .collect();
        let sizes = labels.iter().map(|l| size(l)).collect();
        let mut seqs = Vec::with_capacity(labels.len());
        let mut at = 0;
        for &length in &udb.lengths {
            seqs.push(udb.sequences[at..at + length as usize].to_vec());
            at += length as usize;
        }
        let index = Index::from_parts(udb.header.word_len(), &udb.counts, list);
        return Database::new(labels, seqs, sizes, index);
    }

    let mut reader = Reader::new(fin);
    let (mut labels, mut seqs, mut sizes) = (Vec::new(), Vec::new(), Vec::new());
    let mut record = Record::new();
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            break;
        }
        let header = record.header();
        let abundance = checked_abundance(header)
            .unwrap_or_else(|e| fail(&format!("{e} in read {header} of the database")));
        labels.push(header.to_string());
        seqs.push(record.seq().to_vec());
        sizes.push(abundance as i64);
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    let index = Index::new(word_len, seqs.iter().map(Vec::as_slice));
    if index.without_words() > 0 {
        warn(
            format!(
                "{} of {} sequences give no word to index: they are shorter than the word \
                 length ({word_len}), or each of their words holds a lower-case residue or a \
                 symbol other than A, C, G, T or U",
                index.without_words(),
                seqs.len()
            )
            .as_bytes(),
        );
    }
    Database::new(labels, seqs, sizes, index)
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
    let matches = search_args(command)
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
    };

    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let db = database(matches.get_one::<String>("db").unwrap(), word_len);
    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    let mut searcher = Searcher::new(&db, criteria, scoring, limits);
    let mut minus = Vec::new();
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            break;
        }
        let label = record.header();
        let qsize = checked_abundance(label)
            .unwrap_or_else(|e| fail(&format!("{e} in read {label}"))) as i64;
        let query = record.seq();
        let mut hits = searcher.search(label, qsize, query, false);
        if both {
            reverse_complement_keeping_case(&mut minus, query);
            hits.extend(searcher.search(label, qsize, &minus, true));
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
        let written = if hits.is_empty() {
            write_no_hit(&mut out, label.as_bytes())
        } else {
            hits.iter().try_for_each(|h| {
                let perfect = h.stats.matches == h.stats.length;
                write_hit(
                    &mut out,
                    h.target as u64,
                    query.len(),
                    h.stats.id[iddef],
                    h.minus,
                    if perfect { b"=" } else { h.cigar.as_bytes() },
                    label.as_bytes(),
                    db.labels[h.target].as_bytes(),
                )
            })
        };
        written.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    out.flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
