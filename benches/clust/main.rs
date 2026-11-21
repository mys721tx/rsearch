//! Benchmarks of clustering.
//!
//! `binary` runs the clust program end to end on 10,000 sequences of 250
//! bases. They are copies of 500 templates with one base in twenty changed,
//! inserted or deleted. It clusters them at --id 0.97, and on both strands
//! at --id 0.9 with --maxaccepts 4. compare/benches/clust times vsearch on
//! the same file, after checking that both write the same centroids and
//! hits.

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
    let clust = env!("CARGO_BIN_EXE_clust");
    let mut group = c.benchmark_group("clust/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let wide = ["--id", "0.9", "--maxaccepts", "4", "--strand", "both"];
    for (name, args) in [("id97", &["--id", "0.97"][..]), ("both", &wide[..])] {
        let rargs = [&["-i", fa, "-o", "/dev/null"][..], args].concat();
        group.bench_function(format!("clust/{name}"), |b| b.iter(|| run(clust, &rargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
