//! Times `vsearch --fastq_mergepairs` and its rsearch pipeline on the files
//! of rsearch's mergpa benchmark, at the thread counts rsearch's benchmark
//! takes (`data::thread_counts`; see `common/vsearch.rs`).

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

criterion_group!(benches, bench_binary);
criterion_main!(benches);
