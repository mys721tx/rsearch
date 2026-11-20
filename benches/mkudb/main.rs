//! Benchmarks of building a UDB database.
//!
//! `binary` runs the mkudb program end to end on a FASTA file of 20,000
//! sequences of 1,400 bases, mutated copies of 50 templates (as a reference
//! database of related genes), at the default word length of 8. The input
//! is upper case, so mkudb's file is that of `vsearch --makeudb_usearch
//! --dbmask none`, which is checked before timing. Set `VSEARCH` to a
//! vsearch executable to time it on the same file, on one thread:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench mkudb
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// Sequences of `LENGTH` nucleotides, copies of `TEMPLATES` templates.
const SEQUENCES: usize = 20_000;
const LENGTH: usize = 1_400;
const TEMPLATES: usize = 50;

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

/// Generates copies of random templates with one base in twenty changed.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut out = String::new();
    for i in 0..SEQUENCES {
        let mut seq = templates[rng.below(TEMPLATES)].clone();
        for base in seq.iter_mut() {
            if rng.below(20) == 0 {
                *base = b"ACGT"[rng.below(4)];
            }
        }
        out.push_str(&format!(
            ">seq{i} ref\n{}\n",
            String::from_utf8(seq).unwrap()
        ));
    }
    out
}

fn tmp(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn vsearch_args<'a>(fa: &'a str, out: &'a str) -> [&'a str; 10] {
    [
        "--makeudb_usearch",
        fa,
        "--output",
        out,
        "--dbmask",
        "none",
        "--notrunclabels",
        "--quiet",
        "--threads",
        "1",
    ]
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let input = tmp("mkudb-bench.fa");
    std::fs::write(&input, generate()).unwrap();
    let fa = input.to_str().unwrap();
    let mkudb = env!("CARGO_BIN_EXE_mkudb");
    let null = "/dev/null";
    let vsearch = std::env::var("VSEARCH").ok();

    // both write the same file
    if let Some(vsearch) = &vsearch {
        let (r, v) = (tmp("mkudb-bench-r.udb"), tmp("mkudb-bench-v.udb"));
        run(mkudb, &["-i", fa, "-o", r.to_str().unwrap()]);
        run(vsearch, &vsearch_args(fa, v.to_str().unwrap()));
        assert!(
            std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
            "mkudb and vsearch wrote different files"
        );
    }

    let mut group = c.benchmark_group("mkudb/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    group.bench_function("mkudb", |b| b.iter(|| run(mkudb, &["-i", fa, "-o", null])));
    if let Some(vsearch) = &vsearch {
        let args = vsearch_args(fa, null);
        translated::bench(&mut group, "vsearch", &args);
        group.bench_function("vsearch", |b| b.iter(|| run(vsearch, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
