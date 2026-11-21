//! Benchmarks of LCA taxonomy.
//!
//! `binary` runs the lcatax program end to end on the hits exsrch finds on
//! both strands for 10,000 reads of 150 bases. The reads are drawn from
//! 20,000 distinct sequences, a quarter of them reverse-complemented and a
//! tenth with one base changed. The database holds those 20,000 sequences
//! and 5,000 duplicates, each with a taxonomy. vsearch writes LCA only
//! beside a global search; compare/benches/lcatax times that, its search
//! included, on the same files.

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
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let hits = tmp("lcatax-bench.uc");
    let hits = hits.to_str().unwrap();
    run(
        env!("CARGO_BIN_EXE_exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let lcatax = env!("CARGO_BIN_EXE_lcatax");
    let mut group = c.benchmark_group("lcatax/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    let rargs = [hits, "--db", db, "-o", "/dev/null", "--cutoff", "0.6"];
    group.bench_function("lcatax", |b| b.iter(|| run(lcatax, &rargs)));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
