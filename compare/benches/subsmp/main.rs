//! Times `vsearch --fastx_subsample` and its rsearch pipelines on the file
//! of rsearch's subsmp benchmark: a tenth of the records, and a tenth of
//! the reads by abundance (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/subsmp/input.rs"]
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

    let mut group = c.benchmark_group(data::group("subsmp/binary"));
    group.throughput(Throughput::Elements(data::count(&input, RECORDS)));
    group.sample_size(10);
    for (name, by_size) in [
        ("records", &[][..]),
        ("by_size", &["--sizein", "--sizeout"]),
    ] {
        let args = [
            &[
                "--fastx_subsample",
                fq,
                "--quiet",
                "--notrunclabels",
                "--minseqlength",
                "0",
                "--fastqout",
                "/dev/null",
                "--fastqout_discarded",
                "/dev/null",
                "--randseed",
                "42",
                "--sample_pct",
                "10",
            ][..],
            by_size,
        ]
        .concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
