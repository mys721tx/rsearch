//! Benchmarks of joining read pairs with a pad.
//!
//! `binary` runs the joinpa program end to end on a file of forward FASTQ
//! reads and a file of their mates. compare/benches/joinpa times vsearch
//! on the same files.

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
    let joinpa = env!("CARGO_BIN_EXE_joinpa");

    let mut group = c.benchmark_group("joinpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    let args = ["-i", fwd, "--reverse", rev, "-o", "/dev/null"];
    group.bench_function("joinpa", |b| b.iter(|| run(joinpa, &args)));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
