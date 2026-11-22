//! Benchmarks of multi-parent chimera detection.
//!
//! `binary` runs the mchime program end to end on 1,000 sequences of
//! 1,000 bases, sorted by decreasing abundance: abundant copies of 50
//! parents and rarer chimeras of two to four of them. compare/benches/mchime
//! times vsearch on the same file, after checking that both write the same
//! tables.
//!
//! `simd` times mchime/denovo at each SIMD level this processor has, chosen
//! with `RSEARCH_SIMD` or listed in `RSEARCH_BENCH_SIMD`. It sets the scalar
//! code the compiler vectorises against the hand-written SSE2, AVX2, or NEON.

use std::path::Path;
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
    let (all, _, _) = input_files(dir);
    let all = all.to_str().unwrap();
    let mchime = env!("CARGO_BIN_EXE_mchime");
    let mut group = c.benchmark_group(data::group("mchime/binary"));
    group.throughput(Throughput::Elements(data::count(all, SEQUENCES)));
    group.sample_size(10);
    let args = ["-i", all, "-o", "/dev/null"];
    group.bench_function("mchime/denovo", |b| b.iter(|| run(mchime, &args)));
    group.finish();
}

/// mchime/denovo at each SIMD level this processor has.
fn bench_simd(c: &mut Criterion) {
    if skip() {
        return;
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (all, _, _) = input_files(dir);
    let all = all.to_str().unwrap();
    let mchime = env!("CARGO_BIN_EXE_mchime");
    let mut group = c.benchmark_group(data::group("mchime/simd"));
    group.throughput(Throughput::Elements(data::count(all, SEQUENCES)));
    group.sample_size(10);
    let args = ["-i", all, "-o", "/dev/null"];
    for level in levels::levels() {
        group.bench_function(format!("mchime/denovo/{}", level.name()), |b| {
            b.iter(|| levels::run_at(level, mchime, &args))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary, bench_simd);
criterion_main!(benches);
