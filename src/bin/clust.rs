use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Reader, Record};
use rsearch::hits::write_hit;
use rsearch::label::{Label, checked_abundance};
use rsearch::options::{search_args, search_options};
use rsearch::search::{Database, Hit, Limits, MIN_WORD_MATCHES, Searcher};
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

/// whether `a` is a better centroid than `b`, as vsearch picks among the
/// accepted hits: the higher identity, then the first target; with
/// `by_size`, the more abundant target first.
fn better(a: &Hit, b: &Hit, iddef: usize, sizes: &[i64], by_size: bool) -> bool {
    if by_size && sizes[a.target] != sizes[b.target] {
        return sizes[a.target] > sizes[b.target];
    }
    let (x, y) = (a.stats.id[iddef], b.stats.id[iddef]);
    if x != y {
        return x > y;
    }
    a.target < b.target
}

fn main() {
    let command = Command::new("clust")
        .about("Cluster sequences around centroids, in input order; write centroids and hits.")
        .long_about(
            "Cluster sequences around centroids, in input order; write centroids and \
             hits.\n\n\
             As vsearch --cluster_smallmem --usersort: each sequence is searched against \
             the centroids so far, as by dbsrch (lower case masked), and joins the cluster \
             of its best accepted hit, or becomes a centroid. The centroids are written \
             with the abundance of their cluster (;size=), in the input format; the hits, \
             with --hits, in UC format: an S record for each centroid and an H record for \
             each member, in input order, then a C record for each cluster. Sorted first \
             by sortby, the input gives --cluster_fast and --cluster_size.",
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
                .help("Output centroids in the input format, default to stdout"),
        )
        .arg(
            Arg::new("hits")
                .long("hits")
                .value_name("HITS")
                .help("Output UC file of the clusters"),
        )
        .arg(
            Arg::new("strand")
                .long("strand")
                .value_name("STRAND")
                .value_parser(["plus", "both"])
                .default_value("plus")
                .help("Search the plus strand of the sequences, or both"),
        )
        .arg(
            Arg::new("maxaccepts")
                .long("maxaccepts")
                .value_name("N")
                .value_parser(value_parser!(u64))
                .default_value("1")
                .help("Stop after this many accepted centroids a strand, 0 for no limit"),
        )
        .arg(
            Arg::new("maxrejects")
                .long("maxrejects")
                .value_name("N")
                .value_parser(value_parser!(u64))
                .default_value("32")
                .help("Stop after this many rejected centroids a strand, 0 for no limit"),
        )
        .arg(
            Arg::new("word-len")
                .long("word-len")
                .value_name("LENGTH")
                .value_parser(
                    value_parser!(u32).range(i64::from(MIN_WORD_LEN)..=i64::from(MAX_WORD_LEN)),
                )
                .default_value("8")
                .help("Length of the words indexed"),
        )
        .arg(
            Arg::new("min-word-matches")
                .long("min-word-matches")
                .value_name("N")
                .value_parser(value_parser!(u32))
                .help("Words a centroid shares with the sequence to be a candidate, default by --word-len"),
        )
        .arg(
            Arg::new("size-order")
                .long("size-order")
                .action(ArgAction::SetTrue)
                .help("Join the most abundant of the centroids accepted, rather than the closest"),
        )
        .arg(
            Arg::new("no-sizein")
                .long("no-sizein")
                .action(ArgAction::SetTrue)
                .help("Count each sequence once, rather than by its abundance (;size=)"),
        )
        .arg(
            Arg::new("no-sizeout")
                .long("no-sizeout")
                .action(ArgAction::SetTrue)
                .help("Do not annotate the centroids with their cluster's abundance"),
        );
    let matches = search_args(command)
        .mut_arg("id", |a| a.required(true))
        .get_matches();
    let (scoring, criteria) = search_options(&matches);
    let iddef = criteria.iddef;
    let both = matches.get_one::<String>("strand").unwrap() == "both";
    let word_len = *matches.get_one::<u32>("word-len").unwrap();
    let limits = Limits {
        max_accepts: *matches.get_one::<u64>("maxaccepts").unwrap() as usize,
        max_rejects: *matches.get_one::<u64>("maxrejects").unwrap() as usize,
        min_word_matches: matches
            .get_one::<u32>("min-word-matches")
            .copied()
            .unwrap_or(MIN_WORD_MATCHES[word_len as usize]),
    };
    let by_size = matches.get_flag("size-order");
    let sizein = !matches.get_flag("no-sizein");

    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => create(f),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut hits_out: Option<Box<dyn Write>> = matches.get_one::<String>("hits").map(|f| create(f));

    // every sequence is read before any is clustered
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
    let mut is_fastq = false;
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
        is_fastq = record.is_fastq();
        labels.push(header.to_string());
        seqs.push(record.seq().to_vec());
        quals.push(record.qual().to_vec());
        sizes.push(size as i64);
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }
    let (_, without_words) = word_counts(word_len, seqs.iter().map(Vec::as_slice));
    if without_words > 0 {
        warn(
            format!(
                "{without_words} of {} sequences give no word to index: they are shorter than \
                 the word length ({word_len}), or each of their words holds a lower-case \
                 residue or a symbol other than A, C, G, T or U",
                seqs.len()
            )
            .as_bytes(),
        );
    }

    let mut db = Database::growing(labels, seqs, sizes, word_len);
    let mut searcher = Searcher::new(&db, criteria, scoring, limits);
    let written = |r: io::Result<()>| r.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    // the cluster of each sequence, and the centroid of each cluster
    let mut cluster_of: Vec<usize> = Vec::with_capacity(db.len());
    let mut centroids: Vec<usize> = Vec::new();
    let mut minus = Vec::new();
    for seqno in 0..db.len() {
        let (label, query, size) = (&db.labels[seqno], &db.seqs[seqno], db.sizes[seqno]);
        let mut hits = searcher.search(&db, label, size, query, false);
        if both {
            reverse_complement_keeping_case(&mut minus, query);
            hits.extend(searcher.search(&db, label, size, &minus, true));
        }
        let mut best: Option<&Hit> = None;
        for h in hits.iter().filter(|h| !h.weak) {
            if best.is_none_or(|b| better(h, b, iddef, &db.sizes, by_size)) {
                best = Some(h);
            }
        }
        if let Some(h) = best {
            let cluster = cluster_of[h.target];
            cluster_of.push(cluster);
            if let Some(out) = &mut hits_out {
                let perfect = h.stats.matches == h.stats.length;
                written(write_hit(
                    out,
                    cluster as u64,
                    query.len(),
                    h.stats.id[iddef],
                    h.minus,
                    if perfect { b"=" } else { h.cigar.as_bytes() },
                    label.as_bytes(),
                    db.labels[h.target].as_bytes(),
                ));
            }
        } else {
            let cluster = centroids.len();
            if let Some(out) = &mut hits_out {
                written(writeln!(
                    out,
                    "S\t{cluster}\t{}\t*\t*\t*\t*\t*\t{label}\t*",
                    query.len()
                ));
            }
            cluster_of.push(cluster);
            centroids.push(seqno);
            db.add(seqno);
        }
    }

    // the abundance of each cluster
    let mut abundance = vec![0u64; centroids.len()];
    for (seqno, &cluster) in cluster_of.iter().enumerate() {
        abundance[cluster] += if sizein { db.sizes[seqno] as u64 } else { 1 };
    }
    if let Some(out) = &mut hits_out {
        for (cluster, &seqno) in centroids.iter().enumerate() {
            written(writeln!(
                out,
                "C\t{cluster}\t{}\t*\t*\t*\t*\t*\t{}\t*",
                abundance[cluster], db.labels[seqno]
            ));
        }
        written(out.flush());
    }
    let label = Label {
        sizeout: !matches.get_flag("no-sizeout"),
        ..Label::default()
    };
    let mut writer = Writer::new(&mut out, !is_fastq, 0);
    for (cluster, &seqno) in centroids.iter().enumerate() {
        let header = label.format(
            &db.labels[seqno],
            &db.seqs[seqno],
            abundance[cluster],
            0,
            None,
        );
        written(writer.write(&header, None, &db.seqs[seqno], &quals[seqno]));
    }
    written(writer.flush());
}
