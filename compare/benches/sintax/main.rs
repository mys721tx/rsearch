//! Times `vsearch --sintax` on both strands and its rsearch pipelines on the
//! files of rsearch's sintax benchmark, at the thread counts rsearch's
//! benchmark takes (`data::thread_counts`). It runs against the FASTA database
//! and its UDB file. It first checks that sintax on the most threads writes the
//! table of vsearch --threads 1 (see `common/vsearch.rs`). vsearch on several
//! threads writes the same lines in the order its threads finish.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/sintax/input.rs"]
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
    let (fa, dbfa, udb) = input_files(&vsearch::tool("mkudb"));
    let fa = fa.to_str().unwrap();
    let sintax = vsearch::tool("sintax");

    let mut group = c.benchmark_group(data::group("sintax/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    for (name, db) in [("fasta", &dbfa), ("udb", &udb)] {
        let db = db.to_str().unwrap();
        let (r, v) = (
            data::output(tmp("sintax-bench-r.tsv")),
            data::output(tmp("sintax-bench-v.tsv")),
        );
        let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
        let vbase = [
            "--sintax",
            fa,
            "--db",
            db,
            "--randseed",
            "1",
            "--strand",
            "both",
            "--quiet",
            "--notrunclabels",
        ];
        let counts = data::thread_counts();
        let most = counts.iter().max().unwrap().to_string();
        // both write the same table
        vsearch::run(
            &sintax,
            &[
                "-i", fa, "--db", db, "-o", r, "--seed", "1", "--strand", "both", "-t", &most,
            ],
        );
        vsearch::run(
            &vsearch,
            &[&vbase[..], &["--threads", "1", "--tabbedout", v]].concat(),
        );
        assert!(
            std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
            "sintax and vsearch wrote different tables"
        );
        for threads in counts {
            let t = threads.to_string();
            let vargs = [&vbase[..], &["--threads", &t, "--tabbedout", v]].concat();
            vsearch::bench(&mut group, format!("vsearch/{name}/t{t}"), &vsearch, &vargs);
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
