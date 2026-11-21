//! Benchmarks of exact search.
//!
//! `binary` runs the exsrch program end to end on 100,000 reads of 150
//! bases, on both strands. The reads are drawn from 20,000 distinct
//! sequences, a quarter of them reverse-complemented and a tenth with one
//! base changed. The database holds those 20,000 sequences and 5,000
//! duplicates. compare/benches/exsrch times vsearch on the same files,
//! after checking that both write the same hits.

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
    if skip() {
        return;
    }
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let exsrch = env!("CARGO_BIN_EXE_exsrch");
    let r = data::output(tmp("exsrch-bench-r.uc"));
    let rargs = [
        "-i",
        fa,
        "--db",
        db,
        "-o",
        r.to_str().unwrap(),
        "--strand",
        "both",
    ];
    let mut group = c.benchmark_group(data::group("exsrch/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    group.bench_function("exsrch", |b| b.iter(|| run(exsrch, &rargs)));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
