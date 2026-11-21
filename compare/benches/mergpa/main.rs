//! Times `vsearch --fastq_mergepairs`, on one thread and on its default
//! number, and its rsearch pipeline on the files of rsearch's mergpa
//! benchmark (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/mergpa/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let (forward, reverse) = input_files(&generate());
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());

    let mut group = c.benchmark_group("mergpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
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
    let one = [&args[..], &["--threads", "1"]].concat();
    vsearch::bench(&mut group, "vsearch/1-thread", &vsearch, &one);
    vsearch::bench(&mut group, "vsearch/threads", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
