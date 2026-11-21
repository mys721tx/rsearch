//! Benchmarks of building a UDB database.
//!
//! `binary` runs the mkudb program end to end on a FASTA file of 20,000
//! sequences of 1,400 bases, mutated copies of 50 templates (as a reference
//! database of related genes), at the default word length of 8. The input
//! is upper case. mkudb's file is therefore that of `vsearch
//! --makeudb_usearch --dbmask none`, which compare/benches/mkudb checks
//! before timing it on the same file.

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
    let fa = input.to_str().unwrap();
    let mkudb = env!("CARGO_BIN_EXE_mkudb");
    let null = "/dev/null";

    let mut group = c.benchmark_group(data::group("mkudb/binary"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
    group.sample_size(10);
    group.bench_function("mkudb", |b| b.iter(|| run(mkudb, &["-i", fa, "-o", null])));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
