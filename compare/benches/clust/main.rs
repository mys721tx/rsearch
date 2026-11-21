//! Times `vsearch --cluster_smallmem --usersort --threads 1` and its
//! rsearch pipelines on the file of rsearch's clust benchmark. It runs at
//! --id 0.97, and on both strands at --id 0.9 with --maxaccepts 4. It first
//! checks that clust writes the same centroids and hits (see
//! `common/vsearch.rs`).

use std::path::{Path, PathBuf};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/clust/input.rs"]
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
    let clust = vsearch::tool("clust");
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group("clust/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let wide = ["--id", "0.9", "--maxaccepts", "4", "--strand", "both"];
    for (name, args) in [("id97", &["--id", "0.97"][..]), ("both", &wide[..])] {
        let vbase = [
            "--cluster_smallmem",
            fa,
            "--usersort",
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
        let path = |name: &str| dir.join(format!("clust-bench-{name}"));
        let (rc, ru, vc, vu) = (path("r.fa"), path("r.uc"), path("v.fa"), path("v.uc"));
        let s = |p: &PathBuf| p.to_str().unwrap().to_string();
        vsearch::run(
            &clust,
            &[&["-i", fa, "-o", &s(&rc), "--hits", &s(&ru)][..], args].concat(),
        );
        vsearch::run(
            &vsearch,
            &[&vbase[..], args, &["--centroids", &s(&vc), "--uc", &s(&vu)]].concat(),
        );
        for (r, v) in [(rc, vc), (ru, vu)] {
            assert!(
                std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
                "clust and vsearch wrote different clusters"
            );
        }
        let vargs = [&vbase[..], args, &["--centroids", "/dev/null"]].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &vargs);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
