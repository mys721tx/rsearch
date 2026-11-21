//! Benchmarks of the read length and quality report.
//!
//! `binary` runs the fqstat program end to end on a file of FASTQ reads of
//! varied lengths. compare/benches/fqstat times vsearch on the same input.

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
    if skip() {
        return;
    }
    let input = input_file();
    let fq = input.to_str().unwrap();
    let fqstat = env!("CARGO_BIN_EXE_fqstat");

    let mut group = c.benchmark_group(data::group("fqstat/binary"));
    group.throughput(Throughput::Elements(data::count(&input, READS)));
    group.sample_size(10);
    group.bench_function("fqstat", |b| {
        b.iter(|| run(fqstat, &["-i", fq, "-o", "/dev/null"]))
    });
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
