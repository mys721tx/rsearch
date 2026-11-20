//! Benchmarks of OTU tables.
//!
//! `binary` runs the otutab program end to end on the hits exsrch finds for
//! 100,000 reads of 150 bases from 20 samples, drawn from 20,000 distinct
//! sequences of 2,000 OTUs, a quarter of them reverse-complemented and a
//! tenth with one base changed, against a database of those 20,000
//! sequences and 5,000 duplicates, on both strands. vsearch has no command
//! that only counts hits, so with `VSEARCH` set to a vsearch executable the
//! closest one is timed alongside: `vsearch --search_exact --threads 1`
//! writing the same table, its exact search included, after checking that
//! both write the same classic table:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench otutab
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const READS: usize = 100_000;
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
            ">d{i};otu=OTU{};tax=d:D{},g:G{}\n{}\n",
            i % 2000,
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
        reads.push_str(&format!(
            ">r{i};sample=S{};size={}\n{}\n",
            i % 20,
            1 + i % 3,
            String::from_utf8_lossy(&seq)
        ));
    }
    let (fa, dbfa) = (tmp("otutab-bench.fa"), tmp("otutab-bench-db.fa"));
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

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let hits = tmp("otutab-bench.uc");
    let hits = hits.to_str().unwrap();
    run(
        env!("CARGO_BIN_EXE_exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let otutab = env!("CARGO_BIN_EXE_otutab");
    let vsearch = std::env::var("VSEARCH").ok();
    let (r, v) = (tmp("otutab-bench-r.tsv"), tmp("otutab-bench-v.tsv"));
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    let vbase = [
        "--search_exact",
        fa,
        "--db",
        db,
        "--strand",
        "both",
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
    ];
    // both write the same classic table
    if let Some(vsearch) = &vsearch {
        run(otutab, &[hits, "--db", db, "-o", r]);
        run(vsearch, &[&vbase[..], &["--otutabout", v]].concat());
        assert!(
            std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
            "otutab and vsearch wrote different tables"
        );
    }
    let mut group = c.benchmark_group("otutab/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, rargs, vargs) in [
        ("classic", vec![], vec!["--otutabout", "/dev/null"]),
        (
            "mothur",
            vec!["--mothur"],
            vec!["--mothur_shared_out", "/dev/null"],
        ),
        ("biom", vec!["--biom"], vec!["--biomout", "/dev/null"]),
    ] {
        let rargs = [&[hits, "--db", db, "-o", "/dev/null"][..], &rargs].concat();
        group.bench_function(format!("otutab/{name}"), |b| b.iter(|| run(otutab, &rargs)));
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

criterion_group!(benches, bench_binary);
criterion_main!(benches);
