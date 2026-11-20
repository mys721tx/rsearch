use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::chimera::{Options, Status, Uchime};
use rsearch::fastx::{Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::search::{Database, read_database};
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

/// a number above `low` (or at least `low`, with `inclusive`).
fn above(low: f64, inclusive: bool) -> impl Fn(&str) -> Result<f64, String> + Clone {
    move |text: &str| match text.parse::<f64>() {
        Ok(x) if x.is_finite() && (x > low || (inclusive && x == low)) => Ok(x),
        _ => Err(format!(
            "expected a number {} {low}",
            if inclusive { "of at least" } else { "above" }
        )),
    }
}

/// A record of the input.
struct Query {
    header: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
    size: i64,
    fastq: bool,
}

fn main() {
    let matches = Command::new("uchime")
        .about("Detect chimeras with UCHIME, de novo or against a database.")
        .long_about(
            "Detect chimeras with UCHIME, de novo or against a database.\n\n\
             As vsearch --uchime_denovo on input sorted by abundance (sortby --size), \
             each sequence is checked against the earlier non-chimeras at least --abskew \
             times as abundant; with --db, as --uchime_ref, against the database. The \
             four parts of a sequence are searched for parents, and the two best explain \
             it as a chimera or not, by UCHIME's score or, with --uchime2, by a model \
             matching it perfectly (--uchime2_denovo; --uchime3_denovo with --abskew 16). \
             The non-chimeras are written to -o, the chimeras and borderline sequences to \
             their files, in the input format, with vsearch's table and alignments.",
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
                .help("Output non-chimeras in the input format, default to stdout"),
        )
        .arg(
            Arg::new("db")
                .long("db")
                .value_name("DATABASE")
                .help("Database of parents, FASTA, FASTQ or UDB, rather than de novo"),
        )
        .arg(
            Arg::new("chimeras")
                .long("chimeras")
                .value_name("FILE")
                .help("Output chimeras"),
        )
        .arg(
            Arg::new("borderline")
                .long("borderline")
                .value_name("FILE")
                .help("Output borderline sequences, scoring at least --minh but not chimeras"),
        )
        .arg(
            Arg::new("table")
                .long("table")
                .value_name("FILE")
                .help("Output vsearch's --uchimeout table"),
        )
        .arg(
            Arg::new("alignments")
                .long("alignments")
                .value_name("FILE")
                .help("Output vsearch's --uchimealns alignments of the chimeras"),
        )
        .arg(
            Arg::new("width")
                .long("width")
                .value_name("WIDTH")
                .value_parser(value_parser!(usize))
                .default_value("80")
                .help("Width of the alignments' rows, 0 for no wrapping"),
        )
        .arg(
            Arg::new("uchime2")
                .long("uchime2")
                .action(ArgAction::SetTrue)
                .conflicts_with_all(["db", "minh", "mindiv", "mindiffs"])
                .help("Call chimeras by UCHIME2: a model matching the sequence perfectly"),
        )
        .arg(
            Arg::new("abskew")
                .long("abskew")
                .value_name("SKEW")
                .value_parser(above(1.0, true))
                .conflicts_with("db")
                .help("De novo, parents at least this many times as abundant, default to 2"),
        )
        .arg(
            Arg::new("minh")
                .long("minh")
                .value_name("SCORE")
                .value_parser(above(0.0, false))
                .help("Minimum score of a chimera, default to 0.28"),
        )
        .arg(
            Arg::new("mindiv")
                .long("mindiv")
                .value_name("PERCENT")
                .value_parser(above(0.0, false))
                .help("Minimum divergence of a chimera from its closest parent, default to 0.8"),
        )
        .arg(
            Arg::new("mindiffs")
                .long("mindiffs")
                .value_name("N")
                .value_parser(value_parser!(i32).range(1..))
                .help("Minimum differences on each side of a chimera's crossover, default to 3"),
        )
        .arg(
            Arg::new("xn")
                .long("xn")
                .value_name("WEIGHT")
                .value_parser(above(1.0, false))
                .default_value("8")
                .help("Weight of the no votes"),
        )
        .arg(
            Arg::new("dn")
                .long("dn")
                .value_name("PSEUDOCOUNT")
                .value_parser(above(0.0, false))
                .default_value("1.4")
                .help("Pseudo-count of the no votes"),
        )
        .get_matches();
    let f = |name: &str| matches.get_one::<f64>(name).copied();
    let db_path = matches.get_one::<String>("db");
    let options = Options {
        uchime2: matches.get_flag("uchime2"),
        abskew: db_path.is_none().then(|| f("abskew").unwrap_or(2.0)),
        xn: f("xn").unwrap(),
        dn: f("dn").unwrap(),
        minh: f("minh").unwrap_or(0.28),
        mindiv: f("mindiv").unwrap_or(0.8),
        mindiffs: matches.get_one::<i32>("mindiffs").copied().unwrap_or(3),
        width: *matches.get_one::<usize>("width").unwrap(),
    };

    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut chimeras = matches.get_one::<String>("chimeras").map(|f| create(f));
    let mut borderline = matches.get_one::<String>("borderline").map(|f| create(f));
    let mut table = matches.get_one::<String>("table").map(|f| create(f));
    let mut alignments = matches.get_one::<String>("alignments").map(|f| create(f));

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    let mut next = |record: &mut Record| -> Option<Query> {
        reader.read(record).unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            return None;
        }
        let header = record.header();
        let size =
            checked_abundance(header).unwrap_or_else(|e| fail(&format!("{e} in read {header}")));
        Some(Query {
            header: header.to_string(),
            seq: record.seq().to_vec(),
            qual: record.qual().to_vec(),
            size: size as i64,
            fastq: record.is_fastq(),
        })
    };

    // de novo, the input is the database, read whole, and grows with the
    // non-chimeras. Against a database, the queries are read in turn
    let (mut db, mut queries) = match db_path {
        Some(path) => (
            read_database(path, 8, &mut |w| warn(w)).unwrap_or_else(|e| fail(&e)),
            None,
        ),
        None => {
            let mut all = Vec::new();
            while let Some(q) = next(&mut record) {
                all.push(q);
            }
            let db = Database::growing(
                all.iter().map(|q| q.header.clone()).collect(),
                all.iter().map(|q| q.seq.clone()).collect(),
                all.iter().map(|q| q.size).collect(),
                8,
            );
            (db, Some(all))
        }
    };
    let mut uchime = Uchime::new(&db, options);
    let mut writers: Option<[Writer<Box<dyn Write>>; 3]> = None;
    let written = |r: io::Result<()>| r.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    let mut seqno = 0;
    loop {
        let query = match &mut queries {
            Some(all) => {
                if seqno == all.len() {
                    break;
                }
                std::mem::replace(
                    &mut all[seqno],
                    Query {
                        header: String::new(),
                        seq: Vec::new(),
                        qual: Vec::new(),
                        size: 0,
                        fastq: false,
                    },
                )
            }
            None => match next(&mut record) {
                Some(q) => q,
                None => break,
            },
        };
        let verdict = uchime.evaluate(&db, &query.header, query.size, &query.seq);
        if let Some(t) = &mut table {
            written(t.write_all(verdict.table.as_bytes()));
        }
        if let (Some(a), Some(text)) = (&mut alignments, &verdict.alignment) {
            written(a.write_all(text.as_bytes()));
        }
        let fastq = query.fastq;
        let ws = writers.get_or_insert_with(|| {
            let sink = || -> Box<dyn Write> { Box::new(io::sink()) };
            [
                Writer::new(std::mem::replace(&mut out, sink()), !fastq, 0),
                Writer::new(chimeras.take().unwrap_or_else(sink), !fastq, 0),
                Writer::new(borderline.take().unwrap_or_else(sink), !fastq, 0),
            ]
        });
        let w = match verdict.status {
            Status::Chimeric => &mut ws[1],
            Status::Suspicious => &mut ws[2],
            _ => &mut ws[0],
        };
        written(w.write(&query.header, None, &query.seq, &query.qual));
        // de novo, the non-chimeras become parents
        if queries.is_some() && verdict.status < Status::Suspicious {
            db.add(seqno);
        }
        seqno += 1;
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    if let Some(ws) = &mut writers {
        for w in ws.iter_mut() {
            written(w.flush());
        }
    }
    for w in [&mut table, &mut alignments].into_iter().flatten() {
        written(w.flush());
    }
    written(out.flush());
}
