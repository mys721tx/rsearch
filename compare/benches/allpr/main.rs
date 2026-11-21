//! Times `vsearch --allpairs_global --threads 1` and its rsearch pipelines
//! on the file of rsearch's allpr benchmark, accepting all and at --id 0.9.
//! It first checks that allpr writes the same hits (see
//! `common/vsearch.rs`).

use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/allpr/input.rs"]
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
    let allpr = vsearch::tool("allpr");
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group("allpr/binary");
    group.throughput(Throughput::Elements(
        (SEQUENCES * (SEQUENCES - 1) / 2) as u64,
    ));
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
            "--threads",
            "1",
        ];
        // both write the same hits
        let (r, v) = (dir.join("allpr-bench-r.uc"), dir.join("allpr-bench-v.uc"));
        vsearch::run(
            &allpr,
            &[&["-i", fa, "-o", r.to_str().unwrap()][..], &rargs].concat(),
        );
        vsearch::run(
            &vsearch,
            &[&vbase[..], &vargs, &["--uc", v.to_str().unwrap()]].concat(),
        );
        assert!(
            std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
            "allpr and vsearch wrote different hits"
        );
        let vargs = [&vbase[..], &vargs, &["--uc", "/dev/null"]].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &vargs);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
