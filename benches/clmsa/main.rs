//! Benchmarks of aligning clusters.
//!
//! `binary` runs the clmsa program end to end on 20,000 sequences of 250
//! bases. They are copies of 200 templates with one base in twenty changed,
//! inserted or deleted. Their clusters at --id 0.9 are made by clust
//! beforehand. It writes the multiple alignments, the consensus sequences
//! and the profiles. compare/benches/clmsa times vsearch against the
//! pipeline of clust then clmsa three times, after checking that both write
//! the same.

use std::path::Path;
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
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let input = input_file(dir);
    let fa = input.to_str().unwrap();
    let (clust, clmsa) = (env!("CARGO_BIN_EXE_clust"), env!("CARGO_BIN_EXE_clmsa"));
    let hits = dir.join("clmsa-bench.uc");
    let hits = hits.to_str().unwrap();
    run(
        clust,
        &["-i", fa, "-o", "/dev/null", "--hits", hits, "--id", "0.9"],
    );
    let formats: [&[&str]; 3] = [&[], &["--consensus"], &["--profile"]];
    let mut group = c.benchmark_group("clmsa/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    for (name, format) in ["msa", "consensus", "profile"].into_iter().zip(formats) {
        let args = [&["--hits", hits, "-i", fa, "-o", "/dev/null"][..], format].concat();
        group.bench_function(format!("clmsa/{name}"), |b| b.iter(|| run(clmsa, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
