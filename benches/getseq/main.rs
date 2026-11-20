//! Benchmarks of selecting records by header.
//!
//! `matches` measures `Selector::matches` on the headers of parsed reads;
//! `binary` runs the getseq program end to end on a FASTQ file: one label,
//! 10,000 labels from a file, a substring, a word, and words of a field from
//! a file. `hits` runs it on the hits exsrch finds for 100,000 reads of 150
//! bases against 25,000 sequences, selecting the queries with hits, those
//! without, and the targets hit. Set `VSEARCH` to a vsearch executable to
//! time `vsearch --fastx_getseq` and `--fastx_getseqs` on the same files,
//! and `vsearch --search_exact --threads 1 --qmask none --dbmask none`
//! writing `--matched`,
//! `--notmatched` and `--dbmatched`, its exact search included, after
//! checking that it selects the same queries:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench getseq
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::select::Selector;

/// Reads of `LENGTH` nucleotides, from `SAMPLES` samples.
const READS: usize = 100_000;
const LENGTH: usize = 250;
const SAMPLES: usize = 100;

/// A small deterministic generator (xorshift64*). Every run benchmarks the
/// same data without extra dependencies.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// The header of read `i`.
fn header(i: usize) -> String {
    format!("read{i};sample=s{};size=1", i % SAMPLES)
}

/// Generates FASTQ reads with sample annotations.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@{}\n{seq}\n+\n{qual}\n", header(i)));
    }
    fastq
}

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_matches(c: &mut Criterion) {
    let reads = parse(&generate());
    let mut group = c.benchmark_group("getseq/matches");
    group.throughput(Throughput::Elements(READS as u64));
    let labels: Vec<String> = (0..READS).step_by(10).map(header).collect();
    for (name, selector) in [
        ("labels", Selector::labels(&labels, false)),
        ("substring", Selector::labels(["SAMPLE=S7;"], true)),
        ("word", Selector::words(["s7"], None)),
        ("field", Selector::words(["s7", "s42"], Some("sample"))),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    black_box(selector.matches(r.header()));
                }
            })
        });
    }
    group.finish();
}

/// The generated input files, written once into the target directory.
struct Inputs {
    fastq: PathBuf,
    labels: PathBuf,
    words: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fastq: dir.join("getseq-bench.fq"),
        labels: dir.join("getseq-bench-labels.txt"),
        words: dir.join("getseq-bench-words.txt"),
    };
    std::fs::write(&inputs.fastq, generate()).unwrap();
    let labels: String = (0..READS).step_by(10).map(|i| header(i) + "\n").collect();
    std::fs::write(&inputs.labels, labels).unwrap();
    let words: String = (0..SAMPLES).step_by(7).map(|s| format!("s{s}\n")).collect();
    std::fs::write(&inputs.words, words).unwrap();
    inputs
}

/// A benchmark case: getseq's arguments, and vsearch's command and
/// arguments.
struct Case<'a> {
    name: &'static str,
    getseq: Vec<&'a str>,
    command: &'static str,
    vsearch: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let labels = inputs.labels.to_str().unwrap();
    let words = inputs.words.to_str().unwrap();
    let label = "read500;sample=s0;size=1";
    vec![
        Case {
            name: "label",
            getseq: vec!["--label", label],
            command: "--fastx_getseq",
            vsearch: vec!["--label", label],
        },
        Case {
            name: "labels",
            getseq: vec!["--labels", labels],
            command: "--fastx_getseqs",
            vsearch: vec!["--labels", labels],
        },
        Case {
            name: "substring",
            getseq: vec!["--label", "SAMPLE=S7;", "--substr"],
            command: "--fastx_getseq",
            vsearch: vec!["--label", "SAMPLE=S7;", "--label_substr_match"],
        },
        Case {
            name: "word",
            getseq: vec!["--word", "s7"],
            command: "--fastx_getseqs",
            vsearch: vec!["--label_word", "s7"],
        },
        Case {
            name: "field_words",
            getseq: vec!["--words", words, "--field", "sample"],
            command: "--fastx_getseqs",
            vsearch: vec!["--label_words", words, "--label_field", "sample"],
        },
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let getseq = env!("CARGO_BIN_EXE_getseq");
    let vsearch = std::env::var("VSEARCH").ok();
    let fq = inputs.fastq.to_str().unwrap();

    let mut group = c.benchmark_group("getseq/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        let args = [&["-i", fq, "-o", "/dev/null"][..], &case.getseq].concat();
        group.bench_function(format!("getseq/{}", case.name), |b| {
            b.iter(|| run(getseq, &args))
        });
        if let Some(vsearch) = &vsearch {
            let args = [
                &[
                    case.command,
                    fq,
                    "--quiet",
                    "--notrunclabels",
                    "--fastqout",
                    "/dev/null",
                ][..],
                &case.vsearch,
            ]
            .concat();
            translated::bench(&mut group, format!("vsearch/{}", case.name), &args);
            group.bench_function(format!("vsearch/{}", case.name), |b| {
                b.iter(|| run(vsearch, &args))
            });
        }
    }
    group.finish();
}

/// writes reads drawn from 20,000 sequences, a tenth of them changed, and
/// a database of those sequences and 5,000 duplicates.
fn hit_files() -> (PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let seqs: Vec<String> = (0..20_000)
        .map(|_| (0..150).map(|_| b"ACGT"[rng.below(4)] as char).collect())
        .collect();
    let mut db = String::new();
    for (i, seq) in seqs.iter().enumerate() {
        db.push_str(&format!(">d{i}\n{seq}\n"));
    }
    for i in 0..5_000 {
        db.push_str(&format!(">dup{i}\n{}\n", seqs[rng.below(seqs.len())]));
    }
    let mut reads = String::new();
    for i in 0..100_000 {
        let mut seq = seqs[rng.below(seqs.len())].clone().into_bytes();
        if rng.below(10) == 0 {
            let at = rng.below(150);
            seq[at] = b"ACGT"[rng.below(4)];
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fa, dbfa) = (dir.join("getseq-hits.fa"), dir.join("getseq-hits-db.fa"));
    std::fs::write(&fa, reads).unwrap();
    std::fs::write(&dbfa, db).unwrap();
    (fa, dbfa)
}

fn bench_hits(c: &mut Criterion) {
    let (fa, db) = hit_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let hits = dir.join("getseq-hits.uc");
    let hits = hits.to_str().unwrap();
    run(
        env!("CARGO_BIN_EXE_exsrch"),
        &["-i", fa, "--db", db, "-o", hits],
    );
    let getseq = env!("CARGO_BIN_EXE_getseq");
    let vsearch = std::env::var("VSEARCH").ok();
    let vbase = [
        "--search_exact",
        fa,
        "--db",
        db,
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
        "--fasta_width",
        "0",
        "--qmask",
        "none",
        "--dbmask",
        "none",
    ];
    // both select the same queries
    if let Some(vsearch) = &vsearch {
        let (r, v) = (dir.join("getseq-hits-r.fa"), dir.join("getseq-hits-v.fa"));
        run(
            getseq,
            &["--hits", hits, "-i", fa, "-o", r.to_str().unwrap()],
        );
        run(
            vsearch,
            &[&vbase[..], &["--matched", v.to_str().unwrap()]].concat(),
        );
        assert!(std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap());
    }
    let mut group = c.benchmark_group("getseq/hits");
    group.throughput(Throughput::Elements(100_000));
    group.sample_size(10);
    for (name, input, args, vargs) in [
        ("matched", fa, vec![], vec!["--matched", "/dev/null"]),
        (
            "notmatched",
            fa,
            vec!["--invert"],
            vec!["--notmatched", "/dev/null"],
        ),
        (
            "dbmatched",
            db,
            vec!["--targets"],
            vec!["--dbmatched", "/dev/null"],
        ),
    ] {
        let args = [&["--hits", hits, "-i", input, "-o", "/dev/null"][..], &args].concat();
        group.bench_function(format!("getseq/{name}"), |b| b.iter(|| run(getseq, &args)));
        if let Some(vsearch) = &vsearch {
            let vargs = [&vbase[..], &vargs].concat();
            translated::bench(&mut group, format!("vsearch/{name}"), &vargs);
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &vargs))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_matches, bench_binary, bench_hits);
criterion_main!(benches);
