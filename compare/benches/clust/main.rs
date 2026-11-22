//! Times `vsearch --cluster_smallmem --usersort` and its rsearch pipelines
//! on the file of rsearch's clust benchmark, at the thread counts rsearch's
//! benchmark takes (`data::thread_counts`). It runs at --id 0.97, and on
//! both strands at --id 0.9 with --maxaccepts 4. It first checks that clust
//! on the most threads writes the centroids and hits of vsearch --threads 1
//! (see `common/vsearch.rs`).

use std::path::{Path, PathBuf};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/clust/input.rs"]
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
    let clust = vsearch::tool("clust");
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group(data::group("clust/binary"));
    group.throughput(Throughput::Elements(data::count(&input, SEQUENCES)));
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
        ];
        let counts = data::thread_counts();
        let most = counts.iter().max().unwrap().to_string();
        // both write the same centroids and hits
        let path = |name: &str| data::output(dir.join(format!("clust-bench-{name}")));
        let (rc, ru, vc, vu) = (path("r.fa"), path("r.uc"), path("v.fa"), path("v.uc"));
        let s = |p: &PathBuf| p.to_str().unwrap().to_string();
        vsearch::run(
            &clust,
            &[
                &["-i", fa, "-o", &s(&rc), "--hits", &s(&ru), "-t", &most][..],
                args,
            ]
            .concat(),
        );
        vsearch::run(
            &vsearch,
            &[
                &vbase[..],
                args,
                &["--threads", "1", "--centroids", &s(&vc), "--uc", &s(&vu)],
            ]
            .concat(),
        );
        for (r, v) in [(rc, vc), (ru, vu)] {
            assert!(
                std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
                "clust and vsearch wrote different clusters"
            );
        }
        for threads in counts {
            let t = threads.to_string();
            let vargs = [
                &vbase[..],
                args,
                &["--threads", &t, "--centroids", "/dev/null"],
            ]
            .concat();
            vsearch::bench(&mut group, format!("vsearch/{name}/t{t}"), &vsearch, &vargs);
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
