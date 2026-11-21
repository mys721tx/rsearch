//! Benchmarks of shuffling.
//!
//! `shuffle` measures `rsearch::random::shuffle` on a deck of a million;
//! `binary` runs the shufl program end to end on files of FASTA and FASTQ
//! reads. compare/benches/shufl times vsearch on the same files.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::random::{Mt19937_64, shuffle};

mod input;
use input::*;

fn bench_shuffle(c: &mut Criterion) {
    let mut group = c.benchmark_group("shufl/shuffle");
    const DECK: usize = 1_000_000;
    group.throughput(Throughput::Elements(DECK as u64));
    let mut deck: Vec<u32> = (0..DECK as u32).collect();
    group.bench_function("mt19937_64", |b| {
        b.iter(|| {
            shuffle(&mut deck, &mut Mt19937_64::new(42));
            black_box(&deck);
        })
    });
    group.finish();
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let shufl = env!("CARGO_BIN_EXE_shufl");
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("shufl/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, file) in [("fasta", fa), ("fastq", fq)] {
        group.bench_function(format!("shufl/{name}"), |b| {
            b.iter(|| run(shufl, &["--seed", "42", "-i", file, "-o", "/dev/null"]))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_shuffle, bench_binary);
criterion_main!(benches);
