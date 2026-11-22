//! Greedy centroid clustering, as vsearch's `cluster.cpp`. Each sequence,
//! in input order, is searched against the centroids so far. It joins the
//! cluster of its best accepted hit, or becomes a centroid whose words are
//! added to the index. Shared by clust and unoise.
//!
//! On several threads a batch of sequences is searched at once against the
//! centroids found before the batch. The batch is then decided in order.
//! A sequence is searched again, against the centroids so far, when a
//! centroid found earlier in the batch could have changed its search (see
//! [`Searcher::would_change`]). Every sequence thus gets the hits it gets on
//! one thread.

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};

use crate::cluster::reverse_complement_keeping_case;
use crate::fastx::{Reader, Record};
use crate::hits::write_hit;
use crate::label::{Label, checked_abundance};
use crate::options::{search_args, search_options};
use crate::search::{Database, Hit, Horizon, Limits, MIN_WORD_MATCHES, Searcher};
use crate::threads::{Workers, threads, threads_args};
use crate::words::{MAX_WORD_LEN, MIN_WORD_LEN, word_counts};
use crate::writer::Writer;

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

/// adds the arguments of the clustering tools: the input, the centroids
/// and the hits, the search's limits, words and criteria, and abundances.
pub fn cluster_args(command: Command) -> Command {
    threads_args(search_args(
        command
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
            ),
    ))
}

/// The searches of a sequence on each strand searched, and how far they
/// went.
type Searched = Vec<(Vec<Hit>, Horizon)>;

/// clusters as the arguments say, writing the centroids and the hits, and
/// exits on an error. With `unoise_alpha`, a hit is accepted only by
/// UNOISE's abundance skew rule.
pub fn run(matches: &ArgMatches, unoise_alpha: Option<f64>) {
    let (scoring, criteria) = search_options(matches);
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
        top_hits: None,
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
    let new_searcher = || {
        let mut searcher = Searcher::new(&db, criteria.clone(), scoring, limits);
        if let Some(alpha) = unoise_alpha {
            searcher.set_unoise(alpha);
        }
        (searcher, Vec::new())
    };
    let threads = threads(matches);
    let workers = Workers::new(threads, new_searcher);
    // a batch of two sequences a thread, one at a time on one thread. A larger
    // batch searches more sequences at once but against older centroids, so
    // that more are searched again. On real amplicons two a thread did best
    let batch = if threads == 1 { 1 } else { 2 * threads };
    let written = |r: io::Result<()>| r.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    // the cluster of each sequence, and the centroid of each cluster
    let mut cluster_of: Vec<usize> = Vec::with_capacity(db.len());
    let mut centroids: Vec<usize> = Vec::new();
    // each strand of a sequence searched, against the centroids so far
    let search =
        |(searcher, minus): &mut (Searcher, Vec<u8>), db: &Database, seqno: usize| -> Searched {
            let (label, query, size) = (&db.labels[seqno], &db.seqs[seqno], db.sizes[seqno]);
            let mut searched = vec![searcher.search_seen(db, label, size, query, false)];
            if both {
                reverse_complement_keeping_case(minus, query);
                searched.push(searcher.search_seen(db, label, size, minus, true));
            }
            searched
        };
    let mut start = 0;
    while start < db.len() {
        let end = (start + batch).min(db.len());
        let batch_searched = {
            let db = &db;
            workers.map((start..end).collect(), |state, seqno| {
                search(state, db, seqno)
            })
        };
        // the centroids found in this batch so far, with their words
        let mut fresh: Vec<(usize, Vec<u32>)> = Vec::new();
        for (seqno, searched) in (start..end).zip(batch_searched) {
            let mut state = workers.state();
            let stale = fresh.iter().any(|(centroid, words)| {
                searched
                    .iter()
                    .any(|(_, horizon)| state.0.would_change(&db, horizon, *centroid, words))
            });
            let searched = if stale {
                search(&mut state, &db, seqno)
            } else {
                searched
            };
            let hits: Vec<Hit> = searched.into_iter().flat_map(|(hits, _)| hits).collect();
            let (label, query) = (&db.labels[seqno], &db.seqs[seqno]);
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
                if end - start > 1 {
                    let words = state.0.distinct_words(&db.seqs[seqno]);
                    fresh.push((seqno, words));
                }
                drop(state);
                db.add(seqno);
            }
        }
        start = end;
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
