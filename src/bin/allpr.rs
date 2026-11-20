use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, ArgMatches, Command, value_parser};

use rsearch::accept::{Criteria, Verdict};
use rsearch::align16::Search16;
use rsearch::alignment::Stats;
use rsearch::fastx::{Reader, Record};
use rsearch::hits::{write_hit, write_no_hit};
use rsearch::label::checked_abundance;
use rsearch::linmem::LinearMemoryAligner;
use rsearch::scoring::Scoring;

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

/// a share from 0 to 1.
fn share(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(x) if (0.0..=1.0).contains(&x) => Ok(x),
        _ => Err("expected a number from 0 to 1".to_string()),
    }
}

/// the arguments of vsearch's acceptance criteria, scores and hit choice,
/// shared by the tools that align.
fn search_args(command: Command) -> Command {
    let number = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("NUMBER")
            .value_parser(value_parser!(f64))
            .help(help)
    };
    let count = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("COUNT")
            .value_parser(value_parser!(i64).range(0..))
            .help(help)
    };
    let flag = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .action(ArgAction::SetTrue)
            .help(help)
    };
    let gaps = |text: &str, open: bool| -> Result<String, String> {
        Scoring::default()
            .set_gaps(text, open)
            .map(|()| text.to_string())
    };
    command
        .arg(
            Arg::new("id")
                .long("id")
                .value_name("SHARE")
                .value_parser(share)
                .help("Accept hits of at least this identity, 0 to 1"),
        )
        .arg(
            Arg::new("iddef")
                .long("iddef")
                .value_name("DEF")
                .value_parser(value_parser!(u8).range(0..=4))
                .default_value("2")
                .help("Identity definition, 0 to 4, as vsearch's"),
        )
        .arg(
            Arg::new("weak-id")
                .long("weak-id")
                .value_name("SHARE")
                .value_parser(share)
                .help("Weak hits: of at least this identity, below --id"),
        )
        .arg(
            Arg::new("maxhits")
                .long("maxhits")
                .value_name("N")
                .value_parser(value_parser!(u64).range(1..))
                .help("Write at most this many hits of a query, default to all"),
        )
        .arg(flag(
            "top-hits-only",
            "Write only the hits of the top identity",
        ))
        .arg(count(
            "max-diffs",
            "Reject hits with more substitutions and interior gap columns",
        ))
        .arg(count("max-subs", "Reject hits with more substitutions"))
        .arg(count("max-gaps", "Reject hits with more interior gaps"))
        .arg(count("min-cols", "Reject hits with fewer interior columns"))
        .arg(flag(
            "left-just",
            "Reject hits with a terminal gap on the left",
        ))
        .arg(flag(
            "right-just",
            "Reject hits with a terminal gap on the right",
        ))
        .arg(number(
            "query-cov",
            "Reject hits covering less of the query, 0 to 1",
        ))
        .arg(number(
            "target-cov",
            "Reject hits covering less of the target, 0 to 1",
        ))
        .arg(number("max-id", "Reject hits of a higher identity, 0 to 1"))
        .arg(number(
            "mid",
            "Reject hits of a lower identity of letter pairs, in percent",
        ))
        .arg(count(
            "id-prefix",
            "Reject targets whose first residues differ from the query's",
        ))
        .arg(count(
            "id-suffix",
            "Reject targets whose last residues differ from the query's",
        ))
        .arg(number(
            "max-qt",
            "Reject targets unless the query's length is at most this times theirs",
        ))
        .arg(number(
            "min-qt",
            "Reject targets unless the query's length is at least this times theirs",
        ))
        .arg(number(
            "min-sl",
            "Reject pairs whose shorter length is below this share of the longer",
        ))
        .arg(number(
            "max-sl",
            "Reject pairs whose shorter length is above this share of the longer",
        ))
        .arg(count(
            "max-qsize",
            "Reject queries of an abundance above this",
        ))
        .arg(count(
            "min-tsize",
            "Reject targets of an abundance below this",
        ))
        .arg(number(
            "min-size-ratio",
            "Reject targets unless the query's abundance is at least this times theirs",
        ))
        .arg(number(
            "max-size-ratio",
            "Reject targets unless the query's abundance is at most this times theirs",
        ))
        .arg(flag("self", "Reject targets with the query's label"))
        .arg(flag("self-id", "Reject targets identical to the query"))
        .arg(
            Arg::new("match")
                .long("match")
                .value_name("SCORE")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i64).range(-32767..=32767))
                .default_value("2")
                .help("Score of a match"),
        )
        .arg(
            Arg::new("mismatch")
                .long("mismatch")
                .value_name("SCORE")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i64).range(-32767..=32767))
                .default_value("-4")
                .help("Score of a mismatch"),
        )
        .arg(
            Arg::new("gapopen")
                .long("gapopen")
                .value_name("PENALTIES")
                .value_parser(move |s: &str| gaps(s, true))
                .help("Gap opening penalties, as vsearch's, default to 20I/2E"),
        )
        .arg(
            Arg::new("gapext")
                .long("gapext")
                .value_name("PENALTIES")
                .value_parser(move |s: &str| gaps(s, false))
                .help("Gap extension penalties, as vsearch's, default to 2I/1E"),
        )
        .arg(flag("n-mismatch", "Count N as a mismatch"))
}

/// the scores and criteria the arguments give.
fn search_options(matches: &ArgMatches) -> (Scoring, Criteria) {
    let mut scoring = Scoring {
        match_score: *matches.get_one::<i64>("match").unwrap(),
        mismatch: *matches.get_one::<i64>("mismatch").unwrap(),
        n_mismatch: matches.get_flag("n-mismatch"),
        ..Scoring::default()
    };
    if let Some(text) = matches.get_one::<String>("gapopen") {
        scoring.set_gaps(text, true).unwrap();
    }
    if let Some(text) = matches.get_one::<String>("gapext") {
        scoring.set_gaps(text, false).unwrap();
    }
    let mut c = Criteria::default();
    let f = |name: &str| matches.get_one::<f64>(name).copied();
    let i = |name: &str| matches.get_one::<i64>(name).copied();
    if let Some(v) = f("id") {
        c.id = v;
    }
    if let Some(v) = f("weak-id") {
        c.weak_id = v;
    }
    // as vsearch's, weak hits are never above --id
    if c.weak_id > c.id {
        c.weak_id = c.id;
    }
    c.iddef = usize::from(*matches.get_one::<u8>("iddef").unwrap());
    c.max_diffs = i("max-diffs").unwrap_or(c.max_diffs);
    c.max_subs = i("max-subs").unwrap_or(c.max_subs);
    c.max_gaps = i("max-gaps").unwrap_or(c.max_gaps);
    c.min_cols = i("min-cols").unwrap_or(c.min_cols);
    c.left_just = matches.get_flag("left-just");
    c.right_just = matches.get_flag("right-just");
    c.query_cov = f("query-cov").unwrap_or(c.query_cov);
    c.target_cov = f("target-cov").unwrap_or(c.target_cov);
    c.max_id = f("max-id").unwrap_or(c.max_id);
    c.mid = f("mid").unwrap_or(c.mid);
    c.id_prefix = i("id-prefix").map_or(0, |v| v as usize);
    c.id_suffix = i("id-suffix").map_or(0, |v| v as usize);
    c.max_qt = f("max-qt").unwrap_or(c.max_qt);
    c.min_qt = f("min-qt").unwrap_or(c.min_qt);
    c.min_sl = f("min-sl").unwrap_or(c.min_sl);
    c.max_sl = f("max-sl").unwrap_or(c.max_sl);
    c.sizes.max_query_size = i("max-qsize").unwrap_or(c.sizes.max_query_size);
    c.sizes.min_target_size = i("min-tsize").unwrap_or(c.sizes.min_target_size);
    c.sizes.min_size_ratio = f("min-size-ratio").unwrap_or(c.sizes.min_size_ratio);
    c.sizes.max_size_ratio = f("max-size-ratio").unwrap_or(c.sizes.max_size_ratio);
    c.self_label = matches.get_flag("self");
    c.self_id = matches.get_flag("self-id");
    (scoring, c)
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
