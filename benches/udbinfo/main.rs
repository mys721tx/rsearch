//! Benchmarks of reporting a UDB file's header.
//!
//! `binary` runs the udbinfo program end to end on the UDB file mkudb
//! writes for 20,000 random sequences of 1,400 bases at word length 12, a
//! file of about 170 MB. compare/benches/udbinfo times vsearch on the same
//! file, after checking that both report the same header.

use std::process::{Command, Stdio};

use criterion::{Criterion, criterion_group, criterion_main};

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
    let udb = input_file(env!("CARGO_BIN_EXE_mkudb"));
    let udb = udb.to_str().unwrap();
    let udbinfo = env!("CARGO_BIN_EXE_udbinfo");

    let mut group = c.benchmark_group("udbinfo/binary");
    group.sample_size(20);
    group.bench_function("udbinfo", |b| b.iter(|| run(udbinfo, &["-i", udb])));
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
