//! Times `vsearch --search_exact --threads 1 --qmask none --dbmask none`
//! writing user fields, BLAST tabular and alignments, and its rsearch
//! pipelines, on the files of rsearch's hitfmt benchmark. vsearch has no
//! command that only formats hits. Its exact search is therefore timed too.
//! First checks that hitfmt, on exsrch's hits, writes the same user fields
//! (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/hitfmt/input.rs"]
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
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let hits = data::output(tmp("hitfmt-bench.uc"));
    let hits = hits.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let fields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+raw+caln+qrow";
    let search = [
        "--search_exact",
        fa,
        "--db",
        db,
        "--strand",
        "both",
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
        "--qmask",
        "none",
        "--dbmask",
        "none",
    ];
    // both write the same user fields
    let (r, v) = (
        data::output(tmp("hitfmt-bench-r.tsv")),
        data::output(tmp("hitfmt-bench-v.tsv")),
    );
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    vsearch::run(
        &vsearch::tool("hitfmt"),
        &[hits, "-i", fa, "--db", db, "-o", r, "--userout", fields],
    );
    vsearch::run(
        &vsearch,
        &[&search[..], &["--userout", v, "--userfields", fields]].concat(),
    );
    assert!(
        std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
        "hitfmt and vsearch wrote different user fields"
    );
    let mut group = c.benchmark_group(data::group("hitfmt/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    for (name, output) in [
        (
            "userout",
            vec!["--userout", "/dev/null", "--userfields", fields],
        ),
        ("blast6", vec!["--blast6out", "/dev/null"]),
        ("alignments", vec!["--alnout", "/dev/null"]),
    ] {
        let args = [&search[..], &output].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
