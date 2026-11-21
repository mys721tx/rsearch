//! Benchmarks of scrambling.
//!
//! `scramble` measures `Scrambler::scramble` on parsed reads, with k of 1,
//! 2 and 6; `binary` runs the scrmbl program end to end on a FASTQ file
//! with the same k. compare/benches/scrmbl times vsearch on the same file.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::random::{SplitMix64, substream};
use rsearch::scramble::Scrambler;

mod input;
use input::*;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_scramble(c: &mut Criterion) {
    let reads = parse(&generate());
    let mut group = c.benchmark_group("scrmbl/scramble");
    group.throughput(Throughput::Elements(READS as u64));
    let mut scrambler = Scrambler::default();
    let mut seq = Vec::new();
    for k in [1, 2, 6] {
        group.bench_function(format!("k{k}"), |b| {
            b.iter(|| {
                for (i, r) in reads.iter().enumerate() {
                    seq.clear();
                    seq.extend_from_slice(r.seq());
                    let mut generator = SplitMix64::new(substream(42, i as u64 + 1));
                    scrambler.scramble(&mut seq, k, &mut generator);
                    black_box(&seq);
                }
            })
        });
    }
    group.finish();
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fq = input.to_str().unwrap();
    let scrmbl = env!("CARGO_BIN_EXE_scrmbl");

    let mut group = c.benchmark_group("scrmbl/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for k in ["1", "2", "6"] {
        group.bench_function(format!("scrmbl/k{k}"), |b| {
            b.iter(|| {
                run(
                    scrmbl,
                    &["-i", fq, "-o", "/dev/null", "--seed", "42", "--kmer", k],
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_scramble, bench_binary);
criterion_main!(benches);
