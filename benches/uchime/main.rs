//! Benchmarks of chimera detection.
//!
//! `binary` runs the uchime program end to end de novo on 4,000 sequences
//! of 250 bases, sorted by decreasing abundance. They are abundant variants
//! of 100 parents and rarer two-parent chimeras of them, with one base in a
//! hundred changed. It also runs the same sequences but the parents against
//! a database of the parents. compare/benches/uchime times vsearch on the
//! same files, after checking that both write the same tables.

use std::path::{Path, PathBuf};
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
    let (all, queries, db) = input_files(dir);
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();
    let (all, queries, db) = (s(&all), s(&queries), s(&db));
    let uchime = env!("CARGO_BIN_EXE_uchime");
    let mut group = c.benchmark_group("uchime/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let cases: [(&str, Vec<&str>); 2] = [
        ("denovo", vec!["-i", &all]),
        ("ref", vec!["-i", &queries, "--db", &db]),
    ];
    for (name, rargs) in cases {
        let ra = [&rargs[..], &["-o", "/dev/null"]].concat();
        group.bench_function(format!("uchime/{name}"), |b| b.iter(|| run(uchime, &ra)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
