//! Benchmarks of aligning all pairs.
//!
//! `binary` runs the allpr program end to end on 600 sequences of 250
//! bases. They are copies of 20 templates with one base in twenty changed,
//! inserted or deleted. That makes 179,700 alignments, accepting all and at
//! --id 0.9. compare/benches/allpr times vsearch on the same file, after
//! checking that both write the same hits.

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
    let allpr = env!("CARGO_BIN_EXE_allpr");
    let mut group = c.benchmark_group("allpr/binary");
    group.throughput(Throughput::Elements(
        (SEQUENCES * (SEQUENCES - 1) / 2) as u64,
    ));
    group.sample_size(10);
    for (name, rargs) in [("all", vec!["--accept-all"]), ("id90", vec!["--id", "0.9"])] {
        let rargs = [&["-i", fa, "-o", "/dev/null"][..], &rargs].concat();
        group.bench_function(format!("allpr/{name}"), |b| b.iter(|| run(allpr, &rargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
