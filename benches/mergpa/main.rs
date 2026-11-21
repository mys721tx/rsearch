//! Benchmarks of merging read pairs.
//!
//! `merge` measures `rsearch::merge::Merger` on parsed pairs; `binary` runs
//! the mergpa program end to end on two FASTQ files. The pairs are 250
//! bases read from both ends of fragments of 300 to 450, with Illumina-like
//! scores, so that most overlap by 50 to 200 bases and merge.
//! compare/benches/mergpa times vsearch on the same files.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::merge::{Merger, Settings};

mod input;
use input::*;

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

fn bench_binary(c: &mut Criterion) {
    let (forward, reverse) = input_files(&generate());
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());
    let mergpa = env!("CARGO_BIN_EXE_mergpa");

    let mut group = c.benchmark_group("mergpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    let args = ["-i", fwd, "--reverse", rev, "-o", "/dev/null"];
    group.bench_function("mergpa", |b| b.iter(|| run(mergpa, &args)));
    group.finish();
}

criterion_group!(benches, bench_merge, bench_binary);
criterion_main!(benches);
