//! Benchmarks of LCA taxonomy.
//!
//! `binary` runs the lcatax program end to end on the hits exsrch finds for
//! 10,000 reads of 150 bases, drawn from 20,000 distinct sequences, a
//! quarter of them reverse-complemented and a tenth with one base changed,
//! against a database of those 20,000 sequences and 5,000 duplicates, each
//! with a taxonomy, on both strands. vsearch writes LCA only beside a global
//! search, so with `VSEARCH` set to a vsearch executable `vsearch
//! --usearch_global --id 1 --maxaccepts 0 --maxrejects 0 --threads 1` is
//! timed alongside, its search included, after checking that it finds the
//! same LCA:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench lcatax
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const READS: usize = 10_000;
const DISTINCT: usize = 20_000;
const DUPLICATES: usize = 5_000;
const LENGTH: usize = 150;

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

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn tmp(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

/// writes the reads and the database.
fn input_files() -> (PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let seqs: Vec<Vec<u8>> = (0..DISTINCT)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut db = String::new();
    for (i, seq) in seqs.iter().enumerate() {
        db.push_str(&format!(
            ">d{i};tax=d:D{},p:P{},g:G{}\n{}\n",
            i % 2,
            i % 7,
            i % 2000,
            String::from_utf8_lossy(seq)
        ));
    }
    for i in 0..DUPLICATES {
        let seq = &seqs[rng.below(DISTINCT)];
        db.push_str(&format!(">dup{i}\n{}\n", String::from_utf8_lossy(seq)));
    }
    let mut reads = String::new();
    for i in 0..READS {
        let mut seq = seqs[rng.below(DISTINCT)].clone();
        if rng.below(10) == 0 {
            let at = rng.below(LENGTH);
            seq[at] = b"ACGT"[rng.below(4)];
        }
        if rng.below(4) == 0 {
            seq = reverse_complement(&seq);
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8_lossy(&seq)));
    }
    let (fa, dbfa) = (tmp("lcatax-bench.fa"), tmp("lcatax-bench-db.fa"));
    std::fs::write(&fa, reads).unwrap();
    std::fs::write(&dbfa, db).unwrap();
    (fa, dbfa)
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let hits = tmp("lcatax-bench.uc");
    let hits = hits.to_str().unwrap();
    run(
        env!("CARGO_BIN_EXE_exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let lcatax = env!("CARGO_BIN_EXE_lcatax");
    let vsearch = std::env::var("VSEARCH").ok();
    let (r, v) = (tmp("lcatax-bench-r.tsv"), tmp("lcatax-bench-v.tsv"));
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    let vbase = [
        "--usearch_global",
        fa,
        "--db",
        db,
        "--strand",
        "both",
        "--id",
        "1",
        "--maxaccepts",
        "0",
        "--maxrejects",
        "0",
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
    ];
    // both write the same LCA
    if let Some(vsearch) = &vsearch {
        run(lcatax, &[hits, "--db", db, "-o", r, "--cutoff", "0.6"]);
        run(
            vsearch,
            &[&vbase[..], &["--lcaout", v, "--lca_cutoff", "0.6"]].concat(),
        );
        assert!(
            std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
            "lcatax and vsearch wrote different LCA"
        );
    }
    let mut group = c.benchmark_group("lcatax/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    let rargs = [hits, "--db", db, "-o", "/dev/null", "--cutoff", "0.6"];
    group.bench_function("lcatax", |b| b.iter(|| run(lcatax, &rargs)));
    if let Some(vsearch) = &vsearch {
        let vargs = [
            &vbase[..],
            &["--lcaout", "/dev/null", "--lca_cutoff", "0.6"],
        ]
        .concat();
        group.bench_function("vsearch", |b| b.iter(|| run(vsearch, &vargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
