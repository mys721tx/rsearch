//! Times `vsearch --fastx_uniques --topn 10` and its rsearch pipeline
//! (`derep | topn 10`) on the FASTA file of rsearch's topn benchmark.
//! vsearch has no command of its own for topn (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/topn/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let inputs = input_files();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("topn/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    let args = [
        "--fastx_uniques",
        fa,
        "--quiet",
        "--sizein",
        "--sizeout",
        "--notrunclabels",
        "--minseqlength",
        "0",
        "--fasta_width",
        "0",
        "--topn",
        "10",
        "--fastaout",
        "/dev/null",
    ];
    vsearch::bench(&mut group, "vsearch/fastx_uniques_topn_10", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
