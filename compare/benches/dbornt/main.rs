//! Times `vsearch --orient --threads 1` and its rsearch pipelines on the
//! files of rsearch's dbornt benchmark. It runs against the FASTA database
//! (with `--dbmask soft`, the database as it is) and its UDB file. It first
//! checks that dbornt writes the same table (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/dbornt/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let (fa, dbfa, udb) = input_files(&vsearch::tool("mkudb"));
    let fa = fa.to_str().unwrap();
    let dbornt = vsearch::tool("dbornt");
    let null = "/dev/null";

    let mut group = c.benchmark_group("dbornt/binary");
    group.throughput(Throughput::Elements(READS as u64));
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
        let (r, v) = (tmp("dbornt-bench-r.tsv"), tmp("dbornt-bench-v.tsv"));
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
