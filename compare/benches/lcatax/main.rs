//! Times `vsearch --usearch_global --id 1 --maxaccepts 0 --maxrejects 0
//! --threads 1 --lcaout` and its rsearch pipeline on the files of rsearch's
//! lcatax benchmark. vsearch writes LCA only beside a global search. Its
//! search is therefore timed too. First checks that lcatax, on exsrch's
//! hits, finds the same LCA (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/lcatax/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let hits = tmp("lcatax-bench.uc");
    let hits = hits.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let (r, v) = (tmp("lcatax-bench-r.tsv"), tmp("lcatax-bench-v.tsv"));
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    let vbase = [
        "--usearch_global",
        fa,
        "--db",
        db,
        "--strand",
        "both",
        "--id",
        "1",
        "--maxaccepts",
        "0",
        "--maxrejects",
        "0",
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
    ];
    // both write the same LCA
    vsearch::run(
        &vsearch::tool("lcatax"),
        &[hits, "--db", db, "-o", r, "--cutoff", "0.6"],
    );
    vsearch::run(
        &vsearch,
        &[&vbase[..], &["--lcaout", v, "--lca_cutoff", "0.6"]].concat(),
    );
    assert!(
        std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
        "lcatax and vsearch wrote different LCA"
    );
    let mut group = c.benchmark_group("lcatax/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    let vargs = [
        &vbase[..],
        &["--lcaout", "/dev/null", "--lca_cutoff", "0.6"],
    ]
    .concat();
    vsearch::bench(&mut group, "vsearch", &vsearch, &vargs);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
