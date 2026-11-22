//! Benchmarks of denoising.
//!
//! `binary` runs the unoise program end to end on 10,000 sequences of 250
//! bases, at the defaults and on both strands, at each thread count of
//! `data::thread_counts` (one and one per core, or `RSEARCH_BENCH_THREADS`).
//! They are copies of 300 templates with one base in a hundred changed,
//! inserted or deleted. They come in decreasing abundance (so sorted, as
//! unoise takes them). compare/benches/unoise times vsearch at the same
//! thread counts on the same file, after checking that both write the same
//! centroids and hits.
//!
//! `simd` times unoise/default at each SIMD level this processor has, chosen
//! with `RSEARCH_SIMD` or listed in `RSEARCH_BENCH_SIMD`. It sets the scalar
//! code the compiler vectorises against the hand-written SSE2, AVX2, or NEON.

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
    let unoise = env!("CARGO_BIN_EXE_unoise");
    let mut group = c.benchmark_group(data::group("unoise/binary"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
    group.sample_size(10);
    for (name, args) in [("default", &[][..]), ("both", &["--strand", "both"][..])] {
        for threads in data::thread_counts() {
            let t = threads.to_string();
            let rargs = [&["-i", fa, "-o", "/dev/null", "-t", &t][..], args].concat();
            group.bench_function(format!("unoise/{name}/t{t}"), |b| {
                b.iter(|| run(unoise, &rargs))
            });
        }
    }
    group.finish();
}

/// unoise/default at each SIMD level this processor has, with the defaults on
/// one thread.
fn bench_simd(c: &mut Criterion) {
    if skip() {
        return;
    }
    let input = input_file();
    let fa = input.to_str().unwrap();
    let unoise = env!("CARGO_BIN_EXE_unoise");
    let mut group = c.benchmark_group(data::group("unoise/simd"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
    group.sample_size(10);
    let args = ["-i", fa, "-o", "/dev/null", "-t", "1"];
    for level in levels::levels() {
        group.bench_function(format!("unoise/default/{}", level.name()), |b| {
            b.iter(|| levels::run_at(level, unoise, &args))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary, bench_simd);
criterion_main!(benches);
