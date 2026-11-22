//! Times `vsearch --usearch_global` and its rsearch pipelines on the files of
//! rsearch's dbsrch benchmark, at the thread counts rsearch's benchmark takes
//! (`data::thread_counts`). It runs at --id 0.97, and on both strands with --id
//! 0.9, --maxaccepts 8, and --maxrejects 64. It first checks that dbsrch on the
//! most threads writes the hits of vsearch --threads 1 (see
//! `common/vsearch.rs`). vsearch on several threads writes the same lines in
//! the order its threads finish.

use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/dbsrch/input.rs"]
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
    let (q, db) = input_files();
    let (q, db) = (q.to_str().unwrap(), db.to_str().unwrap());
    let dbsrch = vsearch::tool("dbsrch");
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group(data::group("dbsrch/binary"));
    group.throughput(Throughput::Elements(data::count(q, QUERIES)));
    group.sample_size(10);
    let wide = [
        "--id",
        "0.9",
        "--maxaccepts",
        "8",
        "--maxrejects",
        "64",
        "--strand",
        "both",
    ];
    for (name, args) in [("id97", &["--id", "0.97"][..]), ("both", &wide[..])] {
        let vbase = [
            "--usearch_global",
            q,
            "--db",
            db,
            "--quiet",
            "--notrunclabels",
            "--qmask",
            "soft",
            "--dbmask",
            "soft",
            "--uc_allhits",
        ];
        let counts = data::thread_counts();
        let most = counts.iter().max().unwrap().to_string();
        // both write the same hits
        let (r, v) = (
            data::output(dir.join("dbsrch-bench-r.uc")),
            data::output(dir.join("dbsrch-bench-v.uc")),
        );
        vsearch::run(
            &dbsrch,
            &[
                &["-i", q, "--db", db, "-o", r.to_str().unwrap(), "-t", &most][..],
                args,
            ]
            .concat(),
        );
        vsearch::run(
            &vsearch,
            &[
                &vbase[..],
                args,
                &["--threads", "1", "--uc", v.to_str().unwrap()],
            ]
            .concat(),
        );
        assert!(
            std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
            "dbsrch and vsearch wrote different hits"
        );
        for threads in counts {
            let t = threads.to_string();
            let vargs = [&vbase[..], args, &["--threads", &t, "--uc", "/dev/null"]].concat();
            vsearch::bench(&mut group, format!("vsearch/{name}/t{t}"), &vsearch, &vargs);
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
