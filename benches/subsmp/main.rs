//! Benchmarks of subsampling.
//!
//! `draw` measures `rsearch::random::draw` of half a million reads from a
//! million; `binary` runs the subsmp program end to end on a FASTQ file of
//! dereplicated reads: a tenth of the records, and a tenth of the reads by
//! abundance. compare/benches/subsmp times vsearch on the same file.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::random::{Mt19937_64, draw};

mod input;
use input::*;

fn bench_draw(c: &mut Criterion) {
    let mut group = c.benchmark_group("subsmp/draw");
    let masses = vec![1u64; 1_000_000];
    group.throughput(Throughput::Elements(masses.len() as u64));
    group.bench_function("half_of_a_million", |b| {
        b.iter(|| black_box(draw(&masses, 500_000, &mut Mt19937_64::new(42))))
    });
    group.finish();
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fq = input.to_str().unwrap();
    let subsmp = env!("CARGO_BIN_EXE_subsmp");

    let mut group = c.benchmark_group("subsmp/binary");
    group.throughput(Throughput::Elements(RECORDS as u64));
    group.sample_size(10);
    for (name, by_size) in [("records", &[][..]), ("by_size", &["--by-size"])] {
        let args = [
            &[
                "-i",
                fq,
                "-o",
                "/dev/null",
                "--discarded",
                "/dev/null",
                "--seed",
                "42",
                "--pct",
                "10",
            ][..],
            by_size,
        ]
        .concat();
        group.bench_function(format!("subsmp/{name}"), |b| b.iter(|| run(subsmp, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_draw, bench_binary);
criterion_main!(benches);
