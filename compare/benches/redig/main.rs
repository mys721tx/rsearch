//! Times `vsearch --cut` and its rsearch pipeline on the file of rsearch's
//! redig benchmark, writing the fragments of both strands and the
//! sequences never cut (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/redig/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let input = input_file();
    let fa = input.to_str().unwrap();
    let null = "/dev/null";

    let mut group = c.benchmark_group("redig/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let args = [
        "--cut",
        fa,
        "--cut_pattern",
        "G^AATT_C",
        "--quiet",
        "--fasta_width",
        "0",
        "--fastaout",
        null,
        "--fastaout_rev",
        null,
        "--fastaout_discarded",
        null,
    ];
    vsearch::bench(&mut group, "vsearch", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
