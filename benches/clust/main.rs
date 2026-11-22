//! Benchmarks of clustering.
//!
//! `binary` runs the clust program end to end on 10,000 sequences of 250
//! bases. They are copies of 500 templates with one base in twenty changed,
//! inserted or deleted. It clusters them at --id 0.97, and on both strands
//! at --id 0.9 with --maxaccepts 4, at each thread count of
//! `data::thread_counts` (one and one per core, or `RSEARCH_BENCH_THREADS`).
//! compare/benches/clust times vsearch at the same thread counts on the
//! same file, after checking that both write the same centroids and hits.
//!
//! `simd` times clust/id97 at each SIMD level this processor has, chosen with
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
    let input = input_file();
    let fa = input.to_str().unwrap();
    let clust = env!("CARGO_BIN_EXE_clust");
    let mut group = c.benchmark_group(data::group("clust/binary"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
    group.sample_size(10);
    let wide = ["--id", "0.9", "--maxaccepts", "4", "--strand", "both"];
    for (name, args) in [("id97", &["--id", "0.97"][..]), ("both", &wide[..])] {
        for threads in data::thread_counts() {
            let t = threads.to_string();
            let rargs = [&["-i", fa, "-o", "/dev/null", "-t", &t][..], args].concat();
            group.bench_function(format!("clust/{name}/t{t}"), |b| {
                b.iter(|| run(clust, &rargs))
            });
        }
    }
    group.finish();
}

/// clust/id97 at each SIMD level this processor has, with --id 0.97 on one
/// thread.
fn bench_simd(c: &mut Criterion) {
    if skip() {
        return;
    }
    let input = input_file();
    let fa = input.to_str().unwrap();
    let clust = env!("CARGO_BIN_EXE_clust");
    let mut group = c.benchmark_group(data::group("clust/simd"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
    group.sample_size(10);
    let args = ["-i", fa, "-o", "/dev/null", "-t", "1", "--id", "0.97"];
    for level in levels::levels() {
        group.bench_function(format!("clust/id97/{}", level.name()), |b| {
            b.iter(|| levels::run_at(level, clust, &args))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary, bench_simd);
criterion_main!(benches);
