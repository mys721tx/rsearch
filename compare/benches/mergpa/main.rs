//! Times `vsearch --fastq_mergepairs` and its rsearch pipeline on the files
//! of rsearch's mergpa benchmark, at the thread counts rsearch's benchmark
//! takes (`data::thread_counts`; see `common/vsearch.rs`). `repeats` times
//! both on one thread on the pairs with tandem repeats of rsearch's
//! `repeats` benchmark.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/mergpa/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let (forward, reverse) = pair_files();
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());

    let mut group = c.benchmark_group(data::group("mergpa/binary"));
    group.throughput(Throughput::Elements(data::count(&forward, PAIRS)));
    group.sample_size(10);
    let args = [
        "--fastq_mergepairs",
        fwd,
        "--reverse",
        rev,
        "--quiet",
        "--fastqout",
        "/dev/null",
    ];
    for threads in data::thread_counts() {
        let t = threads.to_string();
        let with = [&args[..], &["--threads", &t]].concat();
        vsearch::bench(&mut group, format!("vsearch/t{t}"), &vsearch, &with);
    }
    group.finish();
}

/// vsearch and its rsearch pipeline on the pairs with tandem repeats of
/// rsearch's `repeats` benchmark, on one thread.
fn bench_repeats(c: &mut Criterion) {
    if skip() {
        return;
    }
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let (forward, reverse) = repeat_files();
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());
    let mut group = c.benchmark_group(data::group("mergpa/repeats"));
    group.throughput(Throughput::Elements((PAIRS / 4) as u64));
    group.sample_size(10);
    let args = [
        "--fastq_mergepairs",
        fwd,
        "--reverse",
        rev,
        "--quiet",
        "--fastqout",
        "/dev/null",
        "--threads",
        "1",
    ];
    vsearch::bench(&mut group, "vsearch/repeats/t1", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary, bench_repeats);
criterion_main!(benches);
