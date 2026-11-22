//! Benchmarks of searching a database.
//!
//! `binary` runs the dbsrch program end to end on 2,000 queries of 250
//! bases against a database of 10,000. They are copies of 200 templates
//! with one base in twenty changed, inserted or deleted. It searches at the
//! defaults with --id 0.97, and on both strands with --id 0.9, --maxaccepts
//! 8 and --maxrejects 64, at each thread count of `data::thread_counts`
//! (one and one per core, or `RSEARCH_BENCH_THREADS`).
//! compare/benches/dbsrch times vsearch at the same thread counts on the
//! same files, after checking that both write the same hits.
//!
//! `simd` times dbsrch/id97 at each SIMD level this processor has, chosen with
//! `RSEARCH_SIMD` or listed in `RSEARCH_BENCH_SIMD`. It sets the scalar code
//! the compiler vectorises against the hand-written SSE2, AVX2, or NEON.

use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
#[path = "../common/levels.rs"]
mod levels;
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
    let (q, db) = input_files();
    let (q, db) = (q.to_str().unwrap(), db.to_str().unwrap());
    let dbsrch = env!("CARGO_BIN_EXE_dbsrch");
    let mut group = c.benchmark_group(data::group("dbsrch/binary"));
    group.throughput(Throughput::Elements(data::count(q, QUERIES)));
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
        for threads in data::thread_counts() {
            let t = threads.to_string();
            let rargs = [
                &["-i", q, "--db", db, "-o", "/dev/null", "-t", &t][..],
                args,
            ]
            .concat();
            group.bench_function(format!("dbsrch/{name}/t{t}"), |b| {
                b.iter(|| run(dbsrch, &rargs))
            });
        }
    }
    group.finish();
}

/// dbsrch/id97 at each SIMD level this processor has, with the defaults and
/// --id 0.97 on one thread.
fn bench_simd(c: &mut Criterion) {
    if skip() {
        return;
    }
    let (q, db) = input_files();
    let (q, db) = (q.to_str().unwrap(), db.to_str().unwrap());
    let dbsrch = env!("CARGO_BIN_EXE_dbsrch");
    let mut group = c.benchmark_group(data::group("dbsrch/simd"));
    group.throughput(Throughput::Elements(data::count(q, QUERIES)));
    group.sample_size(10);
    let args = [
        "-i",
        q,
        "--db",
        db,
        "-o",
        "/dev/null",
        "-t",
        "1",
        "--id",
        "0.97",
    ];
    for level in levels::levels() {
        group.bench_function(format!("dbsrch/id97/{}", level.name()), |b| {
            b.iter(|| levels::run_at(level, dbsrch, &args))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary, bench_simd);
criterion_main!(benches);
