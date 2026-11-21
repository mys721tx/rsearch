//! Times `vsearch --fastx_syncpairs` and its rsearch pipeline on the files
//! of rsearch's syncpa benchmark, writing all four outputs (see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/syncpa/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let (forward, reverse) = input_files();
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());
    let null = "/dev/null";

    let mut group = c.benchmark_group("syncpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    let args = [
        "--fastx_syncpairs",
        fwd,
        "--reverse",
        rev,
        "--quiet",
        "--fastqout",
        null,
        "--fastqout_rev",
        null,
        "--fastqout_orphans",
        null,
        "--fastqout_orphans_rev",
        null,
    ];
    vsearch::bench(&mut group, "vsearch", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
