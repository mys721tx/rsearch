//! Times `vsearch --cluster_unoise --minsize 1 --threads 1` and its rsearch
//! pipelines on the file of rsearch's unoise benchmark, at the defaults and
//! on both strands. It first checks that unoise writes the same centroids
//! and hits (see `common/vsearch.rs`).

use std::path::{Path, PathBuf};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/unoise/input.rs"]
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
    let fa = input.to_str().unwrap();
    let unoise = vsearch::tool("unoise");
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    // vsearch sorts the input by abundance, ties by label, as the
    // translated pipeline does with sortby. The check does the same, which
    // matters for real data, where abundances tie
    let sorted = data::output(dir.join("unoise-bench-sorted.fa"));
    let sorted = sorted.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("sortby"),
        &["--size", "-i", fa, "-o", sorted],
    );
    let mut group = c.benchmark_group(data::group("unoise/binary"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
    group.sample_size(10);
    for (name, args) in [("default", &[][..]), ("both", &["--strand", "both"][..])] {
        let vbase = [
            "--cluster_unoise",
            fa,
            "--minsize",
            "1",
            "--sizein",
            "--sizeout",
            "--quiet",
            "--notrunclabels",
            "--qmask",
            "soft",
            "--fasta_width",
            "0",
            "--threads",
            "1",
        ];
        // both write the same centroids and hits
        let path = |name: &str| data::output(dir.join(format!("unoise-bench-{name}")));
        let (rc, ru, vc, vu) = (path("r.fa"), path("r.uc"), path("v.fa"), path("v.uc"));
        let s = |p: &PathBuf| p.to_str().unwrap().to_string();
        vsearch::run(
            &unoise,
            &[&["-i", sorted, "-o", &s(&rc), "--hits", &s(&ru)][..], args].concat(),
        );
        vsearch::run(
            &vsearch,
            &[&vbase[..], args, &["--centroids", &s(&vc), "--uc", &s(&vu)]].concat(),
        );
        for (r, v) in [(rc, vc), (ru, vu)] {
            assert!(
                std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
                "unoise and vsearch wrote different clusters"
            );
        }
        let vargs = [&vbase[..], args, &["--centroids", "/dev/null"]].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &vargs);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
