//! Benchmarks of denoising.
//!
//! `binary` runs the unoise program end to end on 10,000 sequences of 250
//! bases, at the defaults and on both strands. They are copies of 300
//! templates with one base in a hundred changed, inserted or deleted. They
//! come in decreasing abundance (so sorted, as unoise takes them).
//! compare/benches/unoise times vsearch on the same file, after checking
//! that both write the same centroids and hits.

use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fa = input.to_str().unwrap();
    let unoise = env!("CARGO_BIN_EXE_unoise");
    let mut group = c.benchmark_group("unoise/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    for (name, args) in [("default", &[][..]), ("both", &["--strand", "both"][..])] {
        let rargs = [&["-i", fa, "-o", "/dev/null"][..], args].concat();
        group.bench_function(format!("unoise/{name}"), |b| b.iter(|| run(unoise, &rargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
