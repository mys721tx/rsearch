//! Benchmarks of merging read pairs.
//!
//! `merge` measures `rsearch::merge::Merger` on parsed pairs; `binary` runs
//! the mergpa program end to end on two FASTQ files. The pairs are 250
//! bases read from both ends of fragments of 300 to 450, with Illumina-like
//! scores, so that most overlap by 50 to 200 bases and merge. Set `VSEARCH`
//! to a vsearch executable to time `vsearch --fastq_mergepairs` on the same
//! files for comparison, on one thread and on its default number:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench mergpa
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::merge::{Merger, Settings};

/// Pairs of reads of `LENGTH` nucleotides.
const PAIRS: usize = 100_000;
const LENGTH: usize = 250;

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

/// A read's sequence and scores.
type Read = (Vec<u8>, Vec<u8>);

/// Scores around Q30, from Q2 to Q41, with an error now and then.
fn read(rng: &mut Rng, bases: impl Iterator<Item = u8>) -> Read {
    let mut seq = Vec::with_capacity(LENGTH);
    let mut qual = Vec::with_capacity(LENGTH);
    for base in bases {
        let score = (18 + rng.below(23)).max(2) as u8;
        let base = if rng.below(200) == 0 {
            b"ACGT"[rng.below(4)]
        } else {
            base
        };
        seq.push(base);
        qual.push(33 + score);
    }
    (seq, qual)
}

/// Generates the forward and reverse reads of each fragment.
fn generate() -> Vec<(Read, Read)> {
    let mut rng = Rng(0x5eed);
    (0..PAIRS)
        .map(|_| {
            let length = 300 + rng.below(151);
            let fragment: Vec<u8> = (0..length).map(|_| b"ACGT"[rng.below(4)]).collect();
            let forward = read(&mut rng, fragment[..LENGTH].iter().copied());
            let reverse = read(
                &mut rng,
                fragment[length - LENGTH..].iter().rev().map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                }),
            );
            (forward, reverse)
        })
        .collect()
}

fn fastq(reads: impl Iterator<Item = (usize, Read)>, mate: u8) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, (seq, qual)) in reads {
        out.extend_from_slice(format!("@read{i}/{mate}\n").as_bytes());
        out.extend_from_slice(&seq);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&qual);
        out.push(b'\n');
    }
    out
}

fn input_files(pairs: &[(Read, Read)]) -> (PathBuf, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let paths = (dir.join("mergpa-bench-1.fq"), dir.join("mergpa-bench-2.fq"));
    let forward = fastq(pairs.iter().map(|p| p.0.clone()).enumerate(), 1);
    let reverse = fastq(pairs.iter().map(|p| p.1.clone()).enumerate(), 2);
    std::fs::write(&paths.0, forward).unwrap();
    std::fs::write(&paths.1, reverse).unwrap();
    paths
}

fn settings() -> Settings {
    Settings {
        min_overlap: 10,
        max_diffs: 10,
        max_diff_pct: 100.0,
        stagger: false,
        offset: 33,
        min_score: 0,
        max_score: 93,
        min_score_out: 0,
        max_score_out: 41,
    }
}

fn bench_merge(c: &mut Criterion) {
    let pairs = generate();
    let mut group = c.benchmark_group("mergpa/merge");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    group.bench_function("merge", |b| {
        let mut merger = Merger::new(settings());
        b.iter(|| {
            let mut merged = 0;
            for ((fs, fq), (rs, rq)) in &pairs {
                if merger.merge((fs, fq), (rs, rq)).unwrap().is_none() {
                    merged += 1;
                }
            }
            black_box(merged)
        })
    });
    group.finish();
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let (forward, reverse) = input_files(&generate());
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());
    let mergpa = env!("CARGO_BIN_EXE_mergpa");

    let mut group = c.benchmark_group("mergpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    let args = ["-i", fwd, "--reverse", rev, "-o", "/dev/null"];
    group.bench_function("mergpa", |b| b.iter(|| run(mergpa, &args)));
    if let Ok(vsearch) = std::env::var("VSEARCH") {
        let args = [
            "--fastq_mergepairs",
            fwd,
            "--reverse",
            rev,
            "--quiet",
            "--fastqout",
            "/dev/null",
        ];
        let one = [&args[..], &["--threads", "1"]].concat();
        translated::bench(&mut group, "vsearch/1-thread", &one);
        group.bench_function("vsearch/1-thread", |b| b.iter(|| run(&vsearch, &one)));
        translated::bench(&mut group, "vsearch/threads", &args);
        group.bench_function("vsearch/threads", |b| b.iter(|| run(&vsearch, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_merge, bench_binary);
criterion_main!(benches);
