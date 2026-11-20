//! Benchmarks of exact search.
//!
//! `binary` runs the exsrch program end to end: 100,000 reads of 150 bases
//! drawn from 20,000 distinct sequences, a quarter of them
//! reverse-complemented and a tenth with one base changed, against a
//! database of those 20,000 sequences and 5,000 duplicates, on both
//! strands. Set `VSEARCH` to a vsearch executable to time `vsearch
//! --search_exact --uc_allhits --threads 1` on the same files, after
//! checking that both write the same hits:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench exsrch
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
        db.push_str(&format!(">d{i}\n{}\n", String::from_utf8_lossy(seq)));
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
    let (fa, dbfa) = (tmp("exsrch-bench.fa"), tmp("exsrch-bench-db.fa"));
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
    let exsrch = env!("CARGO_BIN_EXE_exsrch");
    let vsearch = std::env::var("VSEARCH").ok();
    let (r, v) = (tmp("exsrch-bench-r.uc"), tmp("exsrch-bench-v.uc"));
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    let rargs = ["-i", fa, "--db", db, "-o", r, "--strand", "both"];
    let vargs = [
        "--search_exact",
        fa,
        "--db",
        db,
        "--uc",
        v,
        "--uc_allhits",
        "--strand",
        "both",
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
    ];
    // both write the same hits
    if let Some(vsearch) = &vsearch {
        run(exsrch, &rargs);
        run(vsearch, &vargs);
        assert!(
            std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
            "exsrch and vsearch wrote different hits"
        );
    }
    let mut group = c.benchmark_group("exsrch/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    group.bench_function("exsrch", |b| b.iter(|| run(exsrch, &rargs)));
    if let Some(vsearch) = &vsearch {
        translated::bench(&mut group, "vsearch", &vargs);
        group.bench_function("vsearch", |b| b.iter(|| run(vsearch, &vargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
