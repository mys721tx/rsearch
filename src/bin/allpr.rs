use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, Command};

use rsearch::accept::Verdict;
use rsearch::align16::Search16;
use rsearch::alignment::Stats;
use rsearch::fastx::{Reader, Record};
use rsearch::hits::{write_hit, write_no_hit};
use rsearch::label::checked_abundance;
use rsearch::linmem::LinearMemoryAligner;
use rsearch::options::{search_args, search_options};

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

/// A hit, aligned: its target, statistics (with the aligner's score) and
/// CIGAR.
struct Hit {
    target: usize,
    stats: Stats,
    cigar: String,
}

fn main() {
    let command = Command::new("allpr")
        .about("Align every pair of sequences in one file; write hits.")
        .long_about(
            "Align every pair of sequences in one file; write hits.\n\n\
             As vsearch --allpairs_global: each sequence is aligned globally with each \
             one after it, by vsearch's 16-bit aligner or, when a score overflows or \
             the sequences are too long for it, its linear-memory one; the hits \
             accepted by the criteria are written in UC format, by decreasing \
             identity, then target, and an N record for a sequence without hits.",
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
                .help("Output UC file of hits, default to stdout"),
        )
        .arg(
            Arg::new("accept-all")
                .long("accept-all")
                .action(ArgAction::SetTrue)
                .help("Accept every hit, the criteria aside"),
        );
    let matches = search_args(command)
        .group(
            ArgGroup::new("acceptance")
                .args(["id", "accept-all"])
                .required(true),
        )
        .get_matches();
    let (scoring, criteria) = search_options(&matches);
    let accept_all = matches.get_flag("accept-all");
    let max_hits = matches
        .get_one::<u64>("maxhits")
        .map_or(usize::MAX, |&n| n as usize);
    let top_only = matches.get_flag("top-hits-only");

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    let (mut labels, mut seqs, mut sizes) = (Vec::new(), Vec::new(), Vec::new());
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            break;
        }
        let size = checked_abundance(record.header())
            .unwrap_or_else(|e| fail(&format!("{e} in read {}", record.header())));
        labels.push(record.header().to_string());
        seqs.push(record.seq().to_vec());
        sizes.push(size as i64);
    }
    if let Some(warning) = reader.stripped_warning() {
        warn(&warning);
    }

    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut search = Search16::new(&scoring);
    let mut fallback = LinearMemoryAligner::new(&scoring);
    let table = scoring.substitutions();
    for q in 0..seqs.len() {
        let query = &seqs[q];
        let candidates: Vec<usize> = (q + 1..seqs.len())
            .filter(|&t| {
                accept_all
                    || criteria
                        .unaligned(&labels[q], query, sizes[q], &labels[t], &seqs[t], sizes[t])
            })
            .collect();
        let mut hits: Vec<Hit> = Vec::new();
        if !candidates.is_empty() {
            search.set_query(query);
            let targets: Vec<&[u8]> = candidates.iter().map(|&t| seqs[t].as_slice()).collect();
            for (&t, aligned) in candidates.iter().zip(search.align(&targets)) {
                // the 16-bit aligner's counts, or the fallback's statistics
                let (cigar, stats) = if aligned.score == i16::MAX {
                    let cigar = fallback.align(query, &seqs[t]);
                    let stats = Stats::new(&cigar, query, &seqs[t], &scoring, &table);
                    (cigar, stats)
                } else {
                    let stats = Stats::from_counts(
                        &aligned.cigar,
                        query.len(),
                        seqs[t].len(),
                        i64::from(aligned.length),
                        i64::from(aligned.matches),
                        i64::from(aligned.mismatches),
                        i64::from(aligned.gaps),
                        i64::from(aligned.score),
                    );
                    (aligned.cigar, stats)
                };
                let stats = stats.unwrap_or_else(|e| fail(&format!("internal error: {e}")));
                let verdict = if accept_all {
                    Verdict::Accepted
                } else {
                    criteria.aligned(&stats, &cigar, query.len(), seqs[t].len(), &scoring)
                };
                if verdict == Verdict::Accepted {
                    hits.push(Hit {
                        target: t,
                        stats,
                        cigar,
                    });
                }
            }
        }
        // by decreasing identity, then target
        let id = |h: &Hit| h.stats.id[criteria.iddef];
        hits.sort_by(|a, b| id(b).total_cmp(&id(a)).then(a.target.cmp(&b.target)));
        hits.truncate(max_hits);
        if top_only && let Some(best) = hits.first().map(id) {
            hits.retain(|h| id(h) >= best);
        }
        let written = if hits.is_empty() {
            write_no_hit(&mut out, labels[q].as_bytes())
        } else {
            hits.iter().try_for_each(|h| {
                let perfect = h.stats.matches == h.stats.length;
                write_hit(
                    &mut out,
                    h.target as u64,
                    query.len(),
                    id(h),
                    false,
                    if perfect { b"=" } else { h.cigar.as_bytes() },
                    labels[q].as_bytes(),
                    labels[h.target].as_bytes(),
                )
            })
        };
        written.unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    out.flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
