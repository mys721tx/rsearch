use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::fastx::{Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::random;
use rsearch::sintax::{Classifier, Database, query_rng, table_fields};
use rsearch::threads::{Workers, threads, threads_args};
use rsearch::udb::{self, Udb};
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

fn create(path: &str) -> Box<dyn Write> {
    Box::new(BufWriter::new(
        File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))),
    ))
}

/// a share from 0 to 1.
fn share(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(x) if (0.0..=1.0).contains(&x) => Ok(x),
        _ => Err("expected a number from 0 to 1".to_string()),
    }
}

/// reads the database: a UDB file, whose word length wins over
/// `word_len`, or a FASTA or FASTQ file indexed at `word_len`.
fn database(path: &str, word_len: u32) -> Database {
    let file = File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}")));
    let size = file
        .metadata()
        .ok()
        .filter(|m| m.is_file())
        .map(|m| m.len());
    let mut fin = BufReader::new(file);
    let start = fin
        .fill_buf()
        .unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
    if start.starts_with(udb::MAGIC) {
        let mut list = Vec::new();
        let keep = |_: &[u32], _: u64, block: &[u32]| list.extend_from_slice(block);
        let read = match size {
            Some(size) => Udb::read(&mut fin, size, keep),
            None => {
                let mut bytes = Vec::new();
                fin.read_to_end(&mut bytes)
                    .unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
                Udb::read(&mut &bytes[..], bytes.len() as u64, keep)
            }
        };
        let udb = read.unwrap_or_else(|e| fail(&e.to_string()));
        let headers: Vec<Vec<u8>> = udb
            .header_starts
            .windows(2)
            .map(|pair| udb.headers[pair[0] as usize..pair[1] as usize - 1].to_vec())
            .collect();
        // vsearch reads the abundance of every sequence of a database
        for header in &headers {
            if let Err(e) = checked_abundance(&String::from_utf8_lossy(header)) {
                fail(&format!("{e} in the database"));
            }
        }
        let index = Index::from_parts(udb.header.word_len(), &udb.counts, list);
        return Database::new(headers, udb.lengths, index);
    }

    let mut reader = Reader::new(fin);
    let (mut headers, mut seqs) = (Vec::new(), Vec::new());
    let mut record = Record::new();
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            break;
        }
        if let Err(e) = checked_abundance(record.header()) {
            fail(&format!("{e} in read {} of the database", record.header()));
        }
        headers.push(record.header().as_bytes().to_vec());
        seqs.push(record.seq().to_vec());
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
    let lengths = seqs.iter().map(|s| s.len() as u32).collect();
    Database::new(headers, lengths, index)
}

fn main() {
    let matches = threads_args(
        Command::new("sintax")
            .about("Classify sequences by taxonomy with SINTAX.")
            .long_about(
                "Classify sequences by taxonomy with SINTAX.\n\n\
             As vsearch --sintax: a hundred times, 32 words drawn at random from a \
             query pick the database sequence sharing the most of them; with fifty \
             picks or more, each level of the ;tax= annotations takes the name most \
             picks agree on. Writes a table: the query's header, the names with the \
             share of picks agreeing, the strand, and with --cutoff the names whose \
             share reaches it. Lower case residues are masked, in the queries and \
             the database.",
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
                    .help("Output table, default to stdout"),
            )
            .arg(
                Arg::new("db")
                    .long("db")
                    .value_name("FILE")
                    .required(true)
                    .help("Database: FASTA, FASTQ or UDB file with ;tax= annotations"),
            )
            .arg(
                Arg::new("word-len")
                    .long("word-len")
                    .value_name("LEN")
                    .default_value("8")
                    .value_parser(
                        value_parser!(u32).range(i64::from(MIN_WORD_LEN)..=i64::from(MAX_WORD_LEN)),
                    )
                    .help("Length of the words, 3 to 15; a UDB database has its own"),
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
                Arg::new("cutoff")
                    .long("cutoff")
                    .value_name("SHARE")
                    .value_parser(share)
                    .help("Also write the names agreed on by at least this share of picks"),
            )
            .arg(
                Arg::new("random")
                    .long("random")
                    .action(ArgAction::SetTrue)
                    .help(
                        "Break ties at random, rather than for the shortest, then first, sequence",
                    ),
            )
            .arg(
                Arg::new("seed")
                    .long("seed")
                    .value_name("SEED")
                    .allow_negative_numbers(true)
                    .value_parser(value_parser!(i64))
                    .help("Seed of the random numbers, default to one from the system"),
            )
            .arg(
                Arg::new("report")
                    .long("report")
                    .value_name("FILE")
                    .help("Output file of the number of queries classified"),
            ),
    )
    .get_matches();
    let word_len = *matches.get_one::<u32>("word-len").unwrap();
    let both_strands = matches.get_one::<String>("strand").unwrap() == "both";
    let cutoff = matches.get_one::<f64>("cutoff").copied().unwrap_or(0.0);
    let seed = random::seed(matches.get_one::<i64>("seed").copied());

    // the table is created and the database read before the queries are
    // opened, as vsearch does
    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut report = matches.get_one::<String>("report").map(|f| create(f));
    let db = database(matches.get_one::<String>("db").unwrap(), word_len);

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    // a query with its place, or the error that ends the input. Nothing is
    // read after an error, and the queries before it are written first
    let (mut queries, mut classified) = (0u64, 0u64);
    let mut read = 0u64;
    let mut failed = false;
    let next = || -> Option<Result<(u64, String, Vec<u8>), String>> {
        if failed {
            return None;
        }
        let query = match reader.read(&mut record) {
            Err(e) => Err(e.to_string()),
            Ok(()) if record.is_empty() => return None,
            Ok(()) => match checked_abundance(record.header()) {
                Err(e) => Err(format!("{e} in read {}", record.header())),
                Ok(_) => Ok((read, record.header().to_string(), record.seq().to_vec())),
            },
        };
        read += 1;
        failed = query.is_err();
        Some(query)
    };
    let random_ties = matches.get_flag("random");
    let threads = threads(&matches);
    let workers = Workers::new(threads, || Classifier::new(&db, both_strands, random_ties));
    // each query's line, and whether it was classified; its random draws
    // depend only on the seed and its place
    let classify = |classifier: &mut Classifier,
                    query: Result<(u64, String, Vec<u8>), String>|
     -> Result<(Vec<u8>, bool), String> {
        let (place, header, seq) = query?;
        let mut rng = query_rng(seed, place);
        let classification = classifier.classify(&seq, &mut rng);
        let fields = table_fields(&db, &classification, cutoff);
        let mut line = Vec::with_capacity(header.len() + fields.len() + 2);
        line.extend_from_slice(header.as_bytes());
        line.push(b'\t');
        line.extend_from_slice(&fields);
        line.push(b'\n');
        Ok((line, classification.picks.len() >= 50))
    };
    workers.stream(64 * threads, next, classify, |result| match result {
        Ok((line, is_classified)) => {
            if is_classified {
                classified += 1;
            }
            queries += 1;
            out.write_all(&line)
                .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
        }
        Err(message) => {
            // the queries before an error are written
            let _ = out.flush();
            fail(&message);
        }
    });
    out.flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    if let Some(report) = &mut report {
        let mut text = format!("Classified {classified} of {queries} sequences");
        if queries > 0 {
            text.push_str(&format!(
                " ({:.2}%)",
                100.0 * classified as f64 / queries as f64
            ));
        }
        text.push('\n');
        report
            .write_all(text.as_bytes())
            .and_then(|()| report.flush())
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
}
