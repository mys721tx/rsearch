//! Benchmarks of aligning all pairs.
//!
//! `binary` runs the allpr program end to end on 600 sequences of 250
//! bases. They are copies of 20 templates with one base in twenty changed,
//! inserted or deleted. That makes 179,700 alignments, accepting all and at
//! --id 0.9, at each thread count of `data::thread_counts` (one and one per
//! core, or `RSEARCH_BENCH_THREADS`). compare/benches/allpr times vsearch at
//! the same thread counts on the same file, after checking that both write
//! the same hits.

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
    let input = input_file();
    let fa = input.to_str().unwrap();
    let allpr = env!("CARGO_BIN_EXE_allpr");
    let mut group = c.benchmark_group(data::group("allpr/binary"));
    group.throughput(Throughput::Elements({
        let n = data::count(&input, SEQUENCES);
        n * n.saturating_sub(1) / 2
    }));
    group.sample_size(10);
    for (name, rargs) in [("all", vec!["--accept-all"]), ("id90", vec!["--id", "0.9"])] {
        for threads in data::thread_counts() {
            let t = threads.to_string();
            let rargs = [&["-i", fa, "-o", "/dev/null", "-t", &t][..], &rargs].concat();
            group.bench_function(format!("allpr/{name}/t{t}"), |b| {
                b.iter(|| run(allpr, &rargs))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
