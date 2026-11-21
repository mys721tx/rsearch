//! Times `vsearch --scramble` and its rsearch pipelines on the file of
//! rsearch's scrmbl benchmark, with k of 1, 2 and 6 (see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/scrmbl/input.rs"]
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
    let input = input_file();
    let fq = input.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("scrmbl/binary"));
    group.throughput(Throughput::Elements(data::count(&input, READS)));
    group.sample_size(10);
    for k in ["1", "2", "6"] {
        let args = [
            "--scramble",
            fq,
            "--quiet",
            "--notrunclabels",
            "--fastqout",
            "/dev/null",
            "--randseed",
            "42",
            "--scramble_kmer",
            k,
        ];
        vsearch::bench(&mut group, format!("vsearch/k{k}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
