//! Times `vsearch --allpairs_global` and its rsearch pipelines on the file of
//! rsearch's allpr benchmark, accepting all and at --id 0.9, at the thread
//! counts rsearch's benchmark takes (`data::thread_counts`). It first checks
//! that allpr on the most threads writes the hits of vsearch --threads 1 (see
//! `common/vsearch.rs`). vsearch on several threads writes the same lines in
//! the order its threads finish.

use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/allpr/input.rs"]
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
    let allpr = vsearch::tool("allpr");
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group(data::group("allpr/binary"));
    group.throughput(Throughput::Elements({
        let n = data::count(&input, SEQUENCES);
        n * n.saturating_sub(1) / 2
    }));
    group.sample_size(10);
    for (name, rargs, vargs) in [
        ("all", vec!["--accept-all"], vec!["--acceptall"]),
        ("id90", vec!["--id", "0.9"], vec!["--id", "0.9"]),
    ] {
        let vbase = [
            "--allpairs_global",
            fa,
            "--quiet",
            "--notrunclabels",
            "--qmask",
            "none",
        ];
        let counts = data::thread_counts();
        let most = counts.iter().max().unwrap().to_string();
        // both write the same hits
        let (r, v) = (
            data::output(dir.join("allpr-bench-r.uc")),
            data::output(dir.join("allpr-bench-v.uc")),
        );
        vsearch::run(
            &allpr,
            &[
                &["-i", fa, "-o", r.to_str().unwrap(), "-t", &most][..],
                &rargs,
            ]
            .concat(),
        );
        vsearch::run(
            &vsearch,
            &[
                &vbase[..],
                &vargs,
                &["--threads", "1", "--uc", v.to_str().unwrap()],
            ]
            .concat(),
        );
        assert!(
            std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
            "allpr and vsearch wrote different hits"
        );
        for threads in counts {
            let t = threads.to_string();
            let vargs = [&vbase[..], &vargs, &["--threads", &t, "--uc", "/dev/null"]].concat();
            vsearch::bench(&mut group, format!("vsearch/{name}/t{t}"), &vsearch, &vargs);
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
