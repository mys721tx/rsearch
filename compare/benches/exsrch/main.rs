//! Times `vsearch --search_exact --uc_allhits --threads 1` on both strands
//! and its rsearch pipeline on the files of rsearch's exsrch benchmark. It
//! first checks that exsrch writes the same hits (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/exsrch/input.rs"]
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
    let (r, v) = (
        data::output(tmp("exsrch-bench-r.uc")),
        data::output(tmp("exsrch-bench-v.uc")),
    );
    let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
    let vargs = [
        "--search_exact",
        fa,
        "--db",
        db,
        "--uc",
        v,
        "--uc_allhits",
        "--strand",
        "both",
        "--quiet",
        "--notrunclabels",
        "--threads",
        "1",
    ];
    // both write the same hits
    vsearch::run(
        &vsearch::tool("exsrch"),
        &["-i", fa, "--db", db, "-o", r, "--strand", "both"],
    );
    vsearch::run(&vsearch, &vargs);
    assert!(
        std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
        "exsrch and vsearch wrote different hits"
    );
    let mut group = c.benchmark_group(data::group("exsrch/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    vsearch::bench(&mut group, "vsearch", &vsearch, &vargs);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
