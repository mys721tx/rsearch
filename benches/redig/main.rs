//! Benchmarks of cutting sequences at a restriction pattern.
//!
//! `binary` runs the redig program end to end on a FASTA file of random
//! sequences of 1,000 bases. Most hold an EcoRI site (G^AATT_C) planted a
//! few times. It writes the fragments of both strands and the sequences
//! never cut. compare/benches/redig times vsearch on the same file.

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
    let input = input_file();
    let fa = input.to_str().unwrap();
    let redig = env!("CARGO_BIN_EXE_redig");
    let null = "/dev/null";

    let mut group = c.benchmark_group("redig/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let args = [
        "-i",
        fa,
        "--pattern",
        "G^AATT_C",
        "-o",
        null,
        "--rev-out",
        null,
        "--uncut",
        null,
    ];
    group.bench_function("redig", |b| b.iter(|| run(redig, &args)));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
