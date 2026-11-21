//! Benchmarks of splitting clusters.
//!
//! `binary` runs the clspl program end to end on 20,000 sequences of 250
//! bases. They are copies of 200 templates with one base in twenty changed,
//! inserted or deleted. Their clusters at --id 0.97 are made by clust
//! beforehand. compare/benches/clspl times vsearch against the pipeline of
//! clust then clspl, after checking that both write the same files.

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
    let (clust, clspl) = (env!("CARGO_BIN_EXE_clust"), env!("CARGO_BIN_EXE_clspl"));
    let hits = dir.join("clspl-bench.uc");
    let hits = hits.to_str().unwrap();
    run(
        clust,
        &["-i", fa, "-o", "/dev/null", "--hits", hits, "--id", "0.97"],
    );
    let rdir = dir.join("clspl-bench-r");
    let _ = std::fs::remove_dir_all(&rdir);
    std::fs::create_dir_all(&rdir).unwrap();
    let rprefix = rdir.join("c");
    let split = [
        "--hits",
        hits,
        "-i",
        fa,
        "--prefix",
        rprefix.to_str().unwrap(),
    ];
    let mut group = c.benchmark_group("clspl/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    group.bench_function("clspl/split", |b| b.iter(|| run(clspl, &split)));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
