//! Benchmarks of keeping the first records.
//!
//! `binary` runs the topn program end to end on files of FASTQ and FASTA
//! reads, keeping 10 records or all of them. Either way it reads and checks
//! the whole input. vsearch has no command of its own for this;
//! compare/benches/topn times `vsearch --fastx_uniques --topn 10` against
//! its pipeline on the same FASTA file.

use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let topn = env!("CARGO_BIN_EXE_topn");
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("topn/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, file, count) in [
        ("fastq_10", fq, "10"),
        ("fastq_all", fq, "1000000"),
        ("fasta_10", fa, "10"),
        ("fasta_all", fa, "1000000"),
    ] {
        group.bench_function(format!("topn/{name}"), |b| {
            b.iter(|| run(topn, &["-i", file, "-o", "/dev/null", count]))
        });
    }

    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
