//! Benchmarks of keeping one region of every read.
//!
//! `binary` runs the subseq program end to end on files of FASTQ and FASTA
//! reads: a region in the middle, and from a start to the end.
//! compare/benches/subseq times vsearch on the same files.

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
    let subseq = env!("CARGO_BIN_EXE_subseq");
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("subseq/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    let regions: [(&str, &[&str]); 2] = [
        ("middle", &["--start", "50", "--end", "200"]),
        ("to_end", &["--start", "20"]),
    ];
    for (format, file) in [("fastq", fq), ("fasta", fa)] {
        for (region, args) in regions {
            let name = format!("{region}_{format}");
            let args = [&["-i", file, "-o", "/dev/null"][..], args].concat();
            group.bench_function(format!("subseq/{name}"), |b| b.iter(|| run(subseq, &args)));
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
