//! Benchmarks of the per-position quality and expected-error table.
//!
//! `binary` runs the eestat program end to end on a file of FASTQ reads of
//! varied lengths. compare/benches/eestat times vsearch on the same input.

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
    let eestat = env!("CARGO_BIN_EXE_eestat");

    let mut group = c.benchmark_group(data::group("eestat/binary"));
    group.throughput(Throughput::Elements(data::count(&input, READS)));
    group.sample_size(10);
    group.bench_function("eestat", |b| {
        b.iter(|| run(eestat, &["-i", fq, "-o", "/dev/null"]))
    });
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
