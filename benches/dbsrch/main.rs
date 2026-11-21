//! Benchmarks of searching a database.
//!
//! `binary` runs the dbsrch program end to end on 2,000 queries of 250
//! bases against a database of 10,000. They are copies of 200 templates
//! with one base in twenty changed, inserted or deleted. It searches at the
//! defaults with --id 0.97, and on both strands with --id 0.9, --maxaccepts
//! 8 and --maxrejects 64. compare/benches/dbsrch times vsearch on the same
//! files, after checking that both write the same hits.

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
    let (q, db) = input_files();
    let (q, db) = (q.to_str().unwrap(), db.to_str().unwrap());
    let dbsrch = env!("CARGO_BIN_EXE_dbsrch");
    let mut group = c.benchmark_group("dbsrch/binary");
    group.throughput(Throughput::Elements(QUERIES as u64));
    group.sample_size(10);
    let wide = [
        "--id",
        "0.9",
        "--maxaccepts",
        "8",
        "--maxrejects",
        "64",
        "--strand",
        "both",
    ];
    for (name, args) in [("id97", &["--id", "0.97"][..]), ("both", &wide[..])] {
        let rargs = [&["-i", q, "--db", db, "-o", "/dev/null"][..], args].concat();
        group.bench_function(format!("dbsrch/{name}"), |b| b.iter(|| run(dbsrch, &rargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
