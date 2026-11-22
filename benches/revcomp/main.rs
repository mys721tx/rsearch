//! Benchmarks of reverse complementing.
//!
//! `complement` measures `reverse_complement_keeping_case` on parsed
//! reads; `binary` runs the revcomp program end to end on files of FASTQ
//! and FASTA reads. compare/benches/revcomp times vsearch on the same
//! files.
//!
//! `simd` runs revcomp at each SIMD level this processor has, chosen with
//! `RSEARCH_SIMD` or listed in `RSEARCH_BENCH_SIMD`. It sets the scalar code
//! the compiler vectorises against the hand-written AVX2 or NEON.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Reader, Record};

mod input;
#[path = "../common/levels.rs"]
mod levels;
use input::*;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_complement(c: &mut Criterion) {
    if skip() {
        return;
    }
    let reads = parse(&text());
    let mut group = c.benchmark_group(data::group("revcomp/complement"));
    group.throughput(Throughput::Elements(reads.len() as u64));
    let mut dest = Vec::new();
    group.bench_function("keeping_case", |b| {
        b.iter(|| {
            for r in &reads {
                reverse_complement_keeping_case(&mut dest, r.seq());
                black_box(&dest);
            }
        })
    });
    group.finish();
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let inputs = input_files();
    let revcomp = env!("CARGO_BIN_EXE_revcomp");
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("revcomp/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for (name, file) in [("fastq", fq), ("fasta", fa)] {
        group.bench_function(format!("revcomp/{name}"), |b| {
            b.iter(|| run(revcomp, &["-i", file, "-o", "/dev/null"]))
        });
    }
    group.finish();
}

/// revcomp on FASTQ and FASTA at each SIMD level this processor has.
fn bench_simd(c: &mut Criterion) {
    if skip() {
        return;
    }
    let inputs = input_files();
    let revcomp = env!("CARGO_BIN_EXE_revcomp");
    let mut group = c.benchmark_group(data::group("revcomp/simd"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for (name, file) in [("fastq", &inputs.fastq), ("fasta", &inputs.fasta)] {
        let args = ["-i", file.to_str().unwrap(), "-o", "/dev/null"];
        for level in levels::levels() {
            group.bench_function(format!("revcomp/{name}/{}", level.name()), |b| {
                b.iter(|| levels::run_at(level, revcomp, &args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_complement, bench_binary, bench_simd);
criterion_main!(benches);
