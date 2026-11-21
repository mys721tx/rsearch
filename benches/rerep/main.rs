//! Benchmarks of rereplication.
//!
//! `binary` runs the rerep program end to end on files of dereplicated
//! reads, with abundances from 1 to 20: FASTA with and without `;size=1` on
//! the copies, and FASTQ. compare/benches/rerep times vsearch on the FASTA
//! file (vsearch reads only FASTA).

use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

/// A benchmark case: rerep's arguments.
struct Case<'a> {
    name: &'static str,
    rerep: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let null = "/dev/null";
    vec![
        Case {
            name: "fasta",
            rerep: vec!["-i", fa, "-o", null],
        },
        Case {
            name: "fasta_no_sizeout",
            rerep: vec!["-i", fa, "-o", null, "--no-sizeout"],
        },
        Case {
            name: "fastq",
            rerep: vec!["-i", fq, "-o", null],
        },
    ]
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
    let rerep = env!("CARGO_BIN_EXE_rerep");

    let mut group = c.benchmark_group(data::group("rerep/binary"));
    group.throughput(Throughput::Elements(inputs.copies as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("rerep/{}", case.name), |b| {
            b.iter(|| run(rerep, &case.rerep))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
