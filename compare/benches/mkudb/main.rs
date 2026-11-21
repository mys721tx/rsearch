//! Times `vsearch --makeudb_usearch --dbmask none` on one thread and its
//! rsearch pipeline on the file of rsearch's mkudb benchmark, after
//! checking that mkudb writes the same file (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/mkudb/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn vsearch_args<'a>(fa: &'a str, out: &'a str) -> [&'a str; 10] {
    [
        "--makeudb_usearch",
        fa,
        "--output",
        out,
        "--dbmask",
        "none",
        "--notrunclabels",
        "--quiet",
        "--threads",
        "1",
    ]
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let input = input_file();
    let fa = input.to_str().unwrap();

    // both write the same file
    let (r, v) = (
        data::output(tmp("mkudb-bench-r.udb")),
        data::output(tmp("mkudb-bench-v.udb")),
    );
    vsearch::run(
        &vsearch::tool("mkudb"),
        &["-i", fa, "-o", r.to_str().unwrap()],
    );
    vsearch::run(&vsearch, &vsearch_args(fa, v.to_str().unwrap()));
    assert!(
        std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
        "mkudb and vsearch wrote different files"
    );

    let mut group = c.benchmark_group(data::group("mkudb/binary"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
    group.sample_size(10);
    vsearch::bench(
        &mut group,
        "vsearch",
        &vsearch,
        &vsearch_args(fa, "/dev/null"),
    );
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
