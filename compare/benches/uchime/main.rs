//! Times `vsearch --uchime_denovo` and `--uchime_ref` and their rsearch
//! pipelines on the files of rsearch's uchime benchmark, the reference search
//! at the thread counts rsearch's benchmark takes (`data::thread_counts`). It
//! first checks that uchime on the most threads writes the tables of vsearch
//! --threads 1 (see `common/vsearch.rs`). vsearch's --uchime_ref on several
//! threads writes the same lines in the order its threads finish.

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
    ];
    let counts = data::thread_counts();
    let most = counts.iter().max().unwrap().to_string();
    let cases: [(&str, Vec<&str>, Vec<&str>, Vec<usize>); 2] = [
        (
            "denovo",
            vec!["-i", &sorted],
            vec!["--uchime_denovo", &all],
            vec![1],
        ),
        (
            "ref",
            vec!["-i", &queries, "--db", &db, "-t", &most],
            vec!["--uchime_ref", &queries, "--db", &db, "--dbmask", "soft"],
            counts.clone(),
        ),
    ];
    for (name, rargs, vargs, threads) in cases {
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
            &[
                &vargs[..],
                &common,
                &["--threads", "1", "--uchimeout", &s(&v)],
            ]
            .concat(),
        );
        assert!(
            std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
            "uchime and vsearch wrote different tables"
        );
        for threads in threads {
            let t = threads.to_string();
            let va = [
                &vargs[..],
                &common,
                &["--threads", &t, "--nonchimeras", "/dev/null"],
            ]
            .concat();
            let id = if name == "denovo" {
                format!("vsearch/{name}")
            } else {
                format!("vsearch/{name}/t{t}")
            };
            vsearch::bench(&mut group, id, &vsearch, &va);
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
