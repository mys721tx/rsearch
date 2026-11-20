use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, Command, value_parser};

use rsearch::chimera::{LongOptions, MAX_PARENTS, Mchime, Status};
use rsearch::fastx::{Reader, Record};
use rsearch::label::checked_abundance;
use rsearch::search::Database;
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

/// a number from `low` to `high`.
fn within(low: f64, high: f64) -> impl Fn(&str) -> Result<f64, String> + Clone {
    move |text: &str| match text.parse::<f64>() {
        Ok(x) if (low..=high).contains(&x) => Ok(x),
        _ => Err(format!("expected a number from {low} to {high}")),
    }
}

fn main() {
    let matches = Command::new("mchime")
        .about("Detect chimeras of two or more parents in long exact sequences.")
        .long_about(
            "Detect chimeras of two or more parents in long exact sequences.\n\n\
             As vsearch --chimeras_denovo on input sorted by abundance (sortby --size): \
             each sequence is checked against the earlier non-chimeras at least --abskew \
             times as abundant. Its parts are searched for parents, and it is a chimera \
             when the longest segments matching them, taken in turn, cover it whole with \
             two parents or more. The non-chimeras are written to -o and the chimeras to \
             --chimeras, in the input format, with vsearch's table and alignments.",
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
            Arg::new("chimeras")
                .long("chimeras")
                .value_name("FILE")
                .help("Output chimeras"),
        )
        .arg(
            Arg::new("table")
                .long("table")
                .value_name("FILE")
                .help("Output vsearch's --tabbedout table of the chimeras"),
        )
        .arg(
            Arg::new("alignments")
                .long("alignments")
                .value_name("FILE")
                .help("Output vsearch's --alnout alignments of the chimeras"),
        )
        .arg(
            Arg::new("width")
                .long("width")
                .value_name("WIDTH")
                .value_parser(value_parser!(usize))
                .default_value("60")
                .help("Width of the alignments' rows, 0 for no wrapping"),
        )
        .arg(
            Arg::new("parts")
                .long("parts")
                .value_name("N")
                .value_parser(value_parser!(u64).range(2..=100))
                .help("Parts a sequence is searched in, default to one a hundred bases"),
        )
        .arg(
            Arg::new("max-parents")
                .long("max-parents")
                .value_name("N")
                .value_parser(value_parser!(u64).range(2..=MAX_PARENTS as u64))
                .default_value("3")
                .help("Most parents of a chimera"),
        )
        .arg(
            Arg::new("min-segment")
                .long("min-segment")
                .value_name("LENGTH")
                .value_parser(value_parser!(u64).range(1..))
                .default_value("10")
                .help("Shortest segment of a parent"),
        )
        .arg(
            Arg::new("diff-pct")
                .long("diff-pct")
                .value_name("PERCENT")
                .value_parser(within(0.0, 50.0))
                .default_value("0")
                .help("Share of differences a segment may hold, 0 to 50"),
        )
        .arg(
            Arg::new("abskew")
                .long("abskew")
                .value_name("SKEW")
                .value_parser(|text: &str| match text.parse::<f64>() {
                    Ok(x) if x.is_finite() && x >= 1.0 => Ok(x),
                    _ => Err("expected a number of at least 1".to_string()),
                })
                .default_value("1")
                .help("Parents at least this many times as abundant"),
        )
        .get_matches();
    let options = LongOptions {
        parts: matches.get_one::<u64>("parts").map(|&n| n as usize),
        max_parents: *matches.get_one::<u64>("max-parents").unwrap() as usize,
        min_segment: *matches.get_one::<u64>("min-segment").unwrap() as usize,
        diff_pct: *matches.get_one::<f64>("diff-pct").unwrap(),
        abskew: *matches.get_one::<f64>("abskew").unwrap(),
        width: *matches.get_one::<usize>("width").unwrap(),
    };

    let out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let chimeras: Box<dyn Write> = match matches.get_one::<String>("chimeras") {
        Some(f) => create(f),
        None => Box::new(io::sink()),
    };
    let mut table = matches.get_one::<String>("table").map(|f| create(f));
    let mut alignments = matches.get_one::<String>("alignments").map(|f| create(f));

    // every sequence is read before any is checked
    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    let (mut labels, mut seqs, mut quals, mut sizes) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let mut fastq = false;
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            break;
        }
        let header = record.header();
        let size =
            checked_abundance(header).unwrap_or_else(|e| fail(&format!("{e} in read {header}")));
        fastq = record.is_fastq();
        labels.push(header.to_string());
        seqs.push(record.seq().to_vec());
        quals.push(record.qual().to_vec());
        sizes.push(size as i64);
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }

    let mut db = Database::growing(labels, seqs, sizes, 8);
    let mut mchime = Mchime::new(&db, options);
    let mut writers = [
        Writer::new(out, !fastq, 0),
        Writer::new(chimeras, !fastq, 0),
    ];
    let written = |r: io::Result<()>| r.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    for seqno in 0..db.len() {
        let (label, seq, size) = (&db.labels[seqno], &db.seqs[seqno], db.sizes[seqno]);
        let verdict = mchime.evaluate(&db, label, size, seq);
        if let Some(t) = &mut table {
            written(t.write_all(verdict.table.as_bytes()));
        }
        if let (Some(a), Some(text)) = (&mut alignments, &verdict.alignment) {
            written(a.write_all(text.as_bytes()));
        }
        let chimeric = verdict.status == Status::Chimeric;
        written(writers[usize::from(chimeric)].write(label, None, seq, &quals[seqno]));
        // the non-chimeras become parents
        if !chimeric {
            db.add(seqno);
        }
    }
    for w in writers.iter_mut() {
        written(w.flush());
    }
    for w in [&mut table, &mut alignments].into_iter().flatten() {
        written(w.flush());
    }
}
