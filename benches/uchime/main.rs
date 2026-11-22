//! Benchmarks of chimera detection.
//!
//! `binary` runs the uchime program end to end de novo on 4,000 sequences
//! of 250 bases, sorted by decreasing abundance. They are abundant variants
//! of 100 parents and rarer two-parent chimeras of them, with one base in a
//! hundred changed. It also runs the same sequences but the parents against
//! a database of the parents, at each thread count of `data::thread_counts`
//! (one and one per core, or `RSEARCH_BENCH_THREADS`); de novo runs on one
//! thread. compare/benches/uchime times vsearch at the same thread counts
//! on the same files, after checking that both write the same tables.
//!
//! `simd` times uchime/ref at each SIMD level this processor has, chosen with
//! `RSEARCH_SIMD` or listed in `RSEARCH_BENCH_SIMD`. It sets the scalar code
//! the compiler vectorises against the hand-written SSE2, AVX2, or NEON.

use std::path::{Path, PathBuf};
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
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (all, queries, db) = input_files(dir);
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();
    let (all, queries, db) = (s(&all), s(&queries), s(&db));
    let uchime = env!("CARGO_BIN_EXE_uchime");
    let mut group = c.benchmark_group(data::group("uchime/binary"));
    group.throughput(Throughput::Elements(data::count(&all, SEQUENCES)));
    group.sample_size(10);
    let ra = ["-i", &all, "-o", "/dev/null"];
    group.bench_function("uchime/denovo", |b| b.iter(|| run(uchime, &ra)));
    for threads in data::thread_counts() {
        let t = threads.to_string();
        let ra = ["-i", &queries, "--db", &db, "-o", "/dev/null", "-t", &t];
        group.bench_function(format!("uchime/ref/t{t}"), |b| b.iter(|| run(uchime, &ra)));
    }
    group.finish();
}

/// uchime/ref at each SIMD level this processor has, against the database on
/// one thread.
fn bench_simd(c: &mut Criterion) {
    if skip() {
        return;
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (all, queries, db) = input_files(dir);
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();
    let (all, queries, db) = (s(&all), s(&queries), s(&db));
    let uchime = env!("CARGO_BIN_EXE_uchime");
    let mut group = c.benchmark_group(data::group("uchime/simd"));
    group.throughput(Throughput::Elements(data::count(&all, SEQUENCES)));
    group.sample_size(10);
    let args = ["-i", &queries, "--db", &db, "-o", "/dev/null", "-t", "1"];
    for level in levels::levels() {
        group.bench_function(format!("uchime/ref/{}", level.name()), |b| {
            b.iter(|| levels::run_at(level, uchime, &args))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary, bench_simd);
criterion_main!(benches);
