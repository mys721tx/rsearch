//! Benchmarks of the residue and score symbol report.
//!
//! `binary` runs the fqchar program end to end on a file of FASTQ reads,
//! with the default tail length and a long one. compare/benches/fqchar
//! times vsearch on the same input.

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
    let fqchar = env!("CARGO_BIN_EXE_fqchar");

    let mut group = c.benchmark_group(data::group("fqchar/binary"));
    group.throughput(Throughput::Elements(data::count(&input, READS)));
    group.sample_size(10);
    for tail in ["4", "15"] {
        group.bench_function(format!("fqchar/tail_{tail}"), |b| {
            b.iter(|| run(fqchar, &["-i", fq, "-o", "/dev/null", "--tail", tail]))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
