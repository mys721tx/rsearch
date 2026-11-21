//! Times `vsearch --search_exact --threads 1` writing classic, mothur and
//! BIOM tables, and its rsearch pipelines, on the files of rsearch's otutab
//! benchmark. vsearch has no command that only counts hits. Its exact
//! search is therefore timed too. First checks that otutab, on exsrch's
//! hits, writes the same classic table (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/otutab/input.rs"]
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
    let hits = tmp("otutab-bench.uc");
    let hits = hits.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let (r, v) = (tmp("otutab-bench-r.tsv"), tmp("otutab-bench-v.tsv"));
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    let vbase = [
        "--search_exact",
        fa,
        "--db",
        db,
        "--strand",
        "both",
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
    ];
    // both write the same classic table
    vsearch::run(&vsearch::tool("otutab"), &[hits, "--db", db, "-o", r]);
    vsearch::run(&vsearch, &[&vbase[..], &["--otutabout", v]].concat());
    assert!(
        std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
        "otutab and vsearch wrote different tables"
    );
    let mut group = c.benchmark_group("otutab/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, output) in [
        ("classic", vec!["--otutabout", "/dev/null"]),
        ("mothur", vec!["--mothur_shared_out", "/dev/null"]),
        ("biom", vec!["--biomout", "/dev/null"]),
    ] {
        let vargs = [&vbase[..], &output].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &vargs);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
