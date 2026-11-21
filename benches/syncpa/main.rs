//! Benchmarks of pairing the reads of two files by label.
//!
//! `binary` runs the syncpa program end to end on a file of forward FASTQ
//! reads and a file of their mates in shuffled order. A read in ten misses
//! its mate on each side. It writes all four outputs.
//! compare/benches/syncpa times vsearch on the same files.

use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let (forward, reverse) = input_files();
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());
    let syncpa = env!("CARGO_BIN_EXE_syncpa");
    let null = "/dev/null";

    let mut group = c.benchmark_group("syncpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    let args = [
        "-i",
        fwd,
        "--reverse",
        rev,
        "-o",
        null,
        "--rev-out",
        null,
        "--orphans",
        null,
        "--rev-orphans",
        null,
    ];
    group.bench_function("syncpa", |b| b.iter(|| run(syncpa, &args)));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
