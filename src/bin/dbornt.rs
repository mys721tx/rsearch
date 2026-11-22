use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Format, Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::orient::{Orienter, Strand};
use rsearch::threads::{Workers, threads, threads_args};
use rsearch::udb::{self, Header, Udb};
use rsearch::words::{MAX_WORD_LEN, MIN_WORD_LEN, word_counts};
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
    Box::new(BufWriter::new(
        File::create(path).unwrap_or_else(|e| fail(&format!("cannot create {path}: {e}"))),
    ))
}

/// The word counts of the database: read from a UDB file, whose word
/// length wins over `word_len`, or counted from a FASTA or FASTQ file.
fn database(path: &str, word_len: u32) -> (u32, Vec<u32>) {
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
        let mut bytes = Vec::new();
        let size = match size {
            Some(size) => size,
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
        if let Ok(header) = Header::parse(first, size)
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
        let ignore = |_: &[u32], _: u64, _: &[u32]| {};
        let read = if bytes.is_empty() {
            Udb::read(&mut fin, size, ignore)
        } else {
            Udb::read(&mut &bytes[..], size, ignore)
        };
        let udb = read.unwrap_or_else(|e| fail(&e.to_string()));
        // vsearch reads the abundance of every sequence of a database
        for pair in udb.header_starts.windows(2) {
            let header = &udb.headers[pair[0] as usize..pair[1] as usize - 1];
            if let Err(e) = checked_abundance(&String::from_utf8_lossy(header)) {
                fail(&format!("{e} in the database"));
            }
        }
        return (udb.header.word_len(), udb.counts);
    }

    let mut reader = Reader::new(fin);
    let mut seqs = Vec::new();
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
        seqs.push(record.seq().to_vec());
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    let (counts, without_words) = word_counts(word_len, seqs.iter().map(Vec::as_slice));
    if without_words > 0 {
        warn(
            format!(
                "{without_words} of {} sequences give no word to index: they are shorter \
                 than the word length ({word_len}), or each of their words holds a lower-case \
                 residue or a symbol other than A, C, G, T or U",
                seqs.len()
            )
            .as_bytes(),
        );
    }
    (word_len, counts)
}

/// vsearch's summary of the orientations.
fn report(forward: u64, reverse: u64, not_oriented: u64) -> String {
    let queries = forward + reverse + not_oriented;
    let line = |label: &str, count: u64| {
        if queries > 0 {
            format!(
                "{label}{count} ({:.2}%)\n",
                100.0 * count as f64 / queries as f64
            )
        } else {
            format!("{label}{count}\n")
        }
    };
    [
        line("Forward oriented sequences: ", forward),
        line("Reverse oriented sequences: ", reverse),
        line("All oriented sequences:     ", forward + reverse),
        line("Not oriented sequences:     ", not_oriented),
        format!("Total number of sequences:  {queries}\n"),
    ]
    .concat()
}

/// The queries a worker orients at once.
const CHUNK: usize = 1000;

/// A chunk of queries (header, sequence, scores), and the error that ends
/// the input after them, if any.
type Chunk = (Vec<(String, Vec<u8>, Vec<u8>)>, Option<String>);

/// What a worker made of a chunk: the bytes of the oriented records, of those
/// not oriented, and of the table; the queries found forward, reverse, and on
/// neither strand; and the error after them.
struct Oriented {
    oriented: Vec<u8>,
    not_oriented: Vec<u8>,
    table: Vec<u8>,
    counts: [u64; 3],
    error: Option<String>,
}

fn main() {
    let file = |name: &'static str, help: &'static str| {
        Arg::new(name).long(name).value_name("FILE").help(help)
    };
    let matches = threads_args(
        Command::new("dbornt")
            .about("Orient sequences by database words; write them forward.")
            .long_about(
                "Orient sequences by database words; write them forward.\n\n\
             As vsearch --orient: each distinct word of a query votes for a strand when \
             more than eight times as many database sequences hold it on that strand as \
             on the other, and a query is oriented when one strand has a vote and four \
             times the votes of the other. Queries found on the reverse strand are \
             written reverse-complemented. Lower case residues are masked, in the \
             queries and the database. Queries are written in their format, FASTA \
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
                    .help("Output file of the oriented sequences, default to stdout"),
            )
            .arg(
                Arg::new("db")
                    .long("db")
                    .value_name("FILE")
                    .required(true)
                    .help("Database: FASTA, FASTQ or UDB file"),
            )
            .arg(
                Arg::new("word-len")
                    .long("word-len")
                    .value_name("LEN")
                    .default_value("12")
                    .value_parser(
                        value_parser!(u32).range(i64::from(MIN_WORD_LEN)..=i64::from(MAX_WORD_LEN)),
                    )
                    .help("Length of the words, 3 to 15; a UDB database has its own"),
            )
            .arg(file(
                "not-oriented",
                "Output file of the sequences not oriented",
            ))
            .arg(file(
                "tabbedout",
                "Output table: label, strand (+, - or ?) and the votes for each strand",
            ))
            .arg(file("report", "Output file of the summary")),
    )
    .get_matches();
    let word_len = *matches.get_one::<u32>("word-len").unwrap();

    // the queries' format is checked first, then the outputs are created and
    // the database read, as vsearch does
    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let is_fastq = reader.format().unwrap_or_else(|e| fail(&e.to_string())) == Some(Format::Fastq);

    let mut oriented: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut not_oriented = matches.get_one::<String>("not-oriented").map(|f| create(f));
    let mut table = matches.get_one::<String>("tabbedout").map(|f| create(f));
    let mut report_out = matches.get_one::<String>("report").map(|f| create(f));

    let db = matches.get_one::<String>("db").unwrap();
    let (word_len, counts) = database(db, word_len);
    let orienter = Orienter::new(word_len, counts);
    let tabled = table.is_some();
    let sorted_out = not_oriented.is_some();

    let flush_all = |oriented: &mut Box<dyn Write>,
                     not_oriented: &mut Option<Box<dyn Write>>,
                     table: &mut Option<Box<dyn Write>>|
     -> io::Result<()> {
        oriented.flush()?;
        if let Some(w) = not_oriented {
            w.flush()?;
        }
        if let Some(t) = table {
            t.flush()?;
        }
        Ok(())
    };

    // the queries in chunks, up to an error; nothing is read after it
    let mut failed = false;
    let next = || -> Option<Chunk> {
        if failed {
            return None;
        }
        let mut queries = Vec::with_capacity(CHUNK);
        let mut record = Record::new();
        while queries.len() < CHUNK {
            let error = match reader.read(&mut record) {
                Err(e) => Some(e.to_string()),
                Ok(()) if record.is_empty() => {
                    failed = true;
                    break;
                }
                Ok(()) => checked_abundance(record.header())
                    .err()
                    .map(|e| format!("{e} in read {}", record.header())),
            };
            if error.is_some() {
                failed = true;
                return Some((queries, error));
            }
            queries.push((
                record.header().to_string(),
                record.seq().to_vec(),
                record.qual().to_vec(),
            ));
        }
        (!queries.is_empty()).then_some((queries, None))
    };
    // each chunk oriented and written out by a worker, as bytes
    let orient = |orienter: &mut Orienter, (queries, error): Chunk| -> Oriented {
        let (mut oriented, mut not_oriented, mut table) = (Vec::new(), Vec::new(), Vec::new());
        let mut counts = [0u64; 3];
        {
            let mut oriented = Writer::new(&mut oriented, !is_fastq, 0);
            let mut not_oriented = Writer::new(&mut not_oriented, !is_fastq, 0);
            let (mut rc, mut rev_qual) = (Vec::new(), Vec::new());
            for (header, seq, qual) in &queries {
                let (strand, votes_fwd, votes_rev) = orienter.orient(seq);
                match strand {
                    Strand::Forward => {
                        counts[0] += 1;
                        oriented.write(header, None, seq, qual)
                    }
                    Strand::Reverse => {
                        counts[1] += 1;
                        reverse_complement_keeping_case(&mut rc, seq);
                        rev_qual.clear();
                        rev_qual.extend(qual.iter().rev());
                        oriented.write(header, None, &rc, &rev_qual)
                    }
                    Strand::Undetermined => {
                        counts[2] += 1;
                        if sorted_out {
                            not_oriented.write(header, None, seq, qual)
                        } else {
                            Ok(())
                        }
                    }
                }
                .expect("writing to memory");
                if tabled {
                    writeln!(
                        table,
                        "{header}\t{}\t{votes_fwd}\t{votes_rev}",
                        strand.symbol()
                    )
                    .expect("writing to memory");
                }
            }
            oriented.flush().expect("writing to memory");
            not_oriented.flush().expect("writing to memory");
        }
        Oriented {
            oriented,
            not_oriented,
            table,
            counts,
            error,
        }
    };
    let threads = threads(&matches);
    let workers = Workers::new(threads, || orienter.clone());
    let (mut forward, mut reverse, mut undetermined) = (0u64, 0u64, 0u64);
    workers.stream(threads, next, orient, |chunk| {
        let written = oriented
            .write_all(&chunk.oriented)
            .and_then(|()| match &mut not_oriented {
                Some(w) => w.write_all(&chunk.not_oriented),
                None => Ok(()),
            })
            .and_then(|()| match &mut table {
                Some(t) => t.write_all(&chunk.table),
                None => Ok(()),
            });
        written.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
        forward += chunk.counts[0];
        reverse += chunk.counts[1];
        undetermined += chunk.counts[2];
        if let Some(message) = chunk.error {
            // the queries before an error are written
            let _ = flush_all(&mut oriented, &mut not_oriented, &mut table);
            fail(&message);
        }
    });
    flush_all(&mut oriented, &mut not_oriented, &mut table)
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    if let Some(out) = &mut report_out {
        out.write_all(report(forward, reverse, undetermined).as_bytes())
            .and_then(|()| out.flush())
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
}
