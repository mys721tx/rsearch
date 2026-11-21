//! Times `vsearch --fastq_join` and its rsearch pipeline on the files of
//! rsearch's joinpa benchmark (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/joinpa/input.rs"]
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
    let (forward, reverse) = input_files();
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());

    let mut group = c.benchmark_group(data::group("joinpa/binary"));
    group.throughput(Throughput::Elements(data::count(&forward, PAIRS)));
    group.sample_size(10);
    let args = [
        "--fastq_join",
        fwd,
        "--reverse",
        rev,
        "--quiet",
        "--fastqout",
        "/dev/null",
    ];
    vsearch::bench(&mut group, "vsearch", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
