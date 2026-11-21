//! Benchmarks of sorting.
//!
//! `binary` runs the sortby program end to end on files of FASTA and FASTQ
//! reads of varied lengths and abundances, by length and by abundance.
//! compare/benches/sortby times vsearch on the same files.

use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let inputs = input_files();
    let sortby = env!("CARGO_BIN_EXE_sortby");
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("sortby/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for (format, file) in [("fasta", fa), ("fastq", fq)] {
        for key in ["length", "size"] {
            let name = format!("{key}_{format}");
            let flag = format!("--{key}");
            group.bench_function(format!("sortby/{name}"), |b| {
                b.iter(|| run(sortby, &[&flag, "-i", file, "-o", "/dev/null"]))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
