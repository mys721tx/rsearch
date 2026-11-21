//! Times `vsearch --uchime_denovo` and `--uchime_ref` with `--threads 1`
//! and their rsearch pipelines on the files of rsearch's uchime benchmark,
//! after checking that uchime writes the same tables (see
//! `common/vsearch.rs`).

use std::path::{Path, PathBuf};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/uchime/input.rs"]
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
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (all, queries, db) = input_files(dir);
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();
    let (all, queries, db) = (s(&all), s(&queries), s(&db));
    let uchime = vsearch::tool("uchime");
    // vsearch sorts de novo input by abundance, ties by label, as the
    // translated pipeline does with sortby. The check does the same, which
    // matters for real data, where abundances tie
    let sorted = s(&data::output(dir.join("uchime-bench-sorted.fa")));
    vsearch::run(
        &vsearch::tool("sortby"),
        &["--size", "-i", &all, "-o", &sorted],
    );
    let mut group = c.benchmark_group(data::group("uchime/binary"));
    group.throughput(Throughput::Elements(data::count(&all, SEQUENCES)));
    group.sample_size(10);
    let common = [
        "--quiet",
        "--notrunclabels",
        "--qmask",
        "soft",
        "--fasta_width",
        "0",
        "--threads",
        "1",
    ];
    let cases: [(&str, Vec<&str>, Vec<&str>); 2] = [
        ("denovo", vec!["-i", &sorted], vec!["--uchime_denovo", &all]),
        (
            "ref",
            vec!["-i", &queries, "--db", &db],
            vec!["--uchime_ref", &queries, "--db", &db, "--dbmask", "soft"],
        ),
    ];
    for (name, rargs, vargs) in cases {
        // both write the same tables
        let (r, v) = (
            data::output(dir.join("uchime-bench-r.tab")),
            data::output(dir.join("uchime-bench-v.tab")),
        );
        vsearch::run(
            &uchime,
            &[&rargs[..], &["-o", "/dev/null", "--table", &s(&r)]].concat(),
        );
        vsearch::run(
            &vsearch,
            &[&vargs[..], &common, &["--uchimeout", &s(&v)]].concat(),
        );
        assert!(
            std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
            "uchime and vsearch wrote different tables"
        );
        let va = [&vargs[..], &common, &["--nonchimeras", "/dev/null"]].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &va);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
