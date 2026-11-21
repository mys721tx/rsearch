//! Benchmarks of multi-parent chimera detection.
//!
//! `binary` runs the mchime program end to end on 1,000 sequences of
//! 1,000 bases, sorted by decreasing abundance: abundant copies of 50
//! parents and rarer chimeras of two to four of them. compare/benches/mchime
//! times vsearch on the same file, after checking that both write the same
//! tables.

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

criterion_group!(benches, bench_binary);
criterion_main!(benches);
