//! Times `vsearch --orient --threads 1` and its rsearch pipelines on the
//! files of rsearch's dbornt benchmark. It runs against the FASTA database
//! (with `--dbmask soft`, the database as it is) and its UDB file. It first
//! checks that dbornt on the most threads of `data::thread_counts` writes
//! the same table (see `common/vsearch.rs`). vsearch orients on one thread.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/dbornt/input.rs"]
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
    let dbornt = vsearch::tool("dbornt");
    let null = "/dev/null";

    let mut group = c.benchmark_group(data::group("dbornt/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    for (name, db) in [("fasta", &dbfa), ("udb", &udb)] {
        let db = db.to_str().unwrap();
        let mut vargs = vec![
            "--orient",
            fa,
            "--db",
            db,
            "--fastaout",
            null,
            "--quiet",
            "--notrunclabels",
            "--threads",
            "1",
        ];
        if name == "fasta" {
            vargs.extend(["--dbmask", "soft"]);
        }
        // both write the same table
        let (r, v) = (
            data::output(tmp("dbornt-bench-r.tsv")),
            data::output(tmp("dbornt-bench-v.tsv")),
        );
        vsearch::run(
            &dbornt,
            &[
                "-i",
                fa,
                "--db",
                db,
                "-o",
                null,
                "--tabbedout",
                r.to_str().unwrap(),
                "-t",
                &data::thread_counts().iter().max().unwrap().to_string(),
            ],
        );
        vsearch::run(
            &vsearch,
            &[&vargs[..], &["--tabbedout", v.to_str().unwrap()]].concat(),
        );
        assert!(
            std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
            "dbornt and vsearch wrote different tables"
        );
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &vargs);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
