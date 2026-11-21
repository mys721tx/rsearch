//! Times `vsearch --udbinfo` and its rsearch pipeline on the UDB file of
//! rsearch's udbinfo benchmark, after checking that udbinfo reports the
//! header vsearch logs (see `common/vsearch.rs`).

use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};

#[path = "../../../benches/udbinfo/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let udb = input_file(&vsearch::tool("mkudb"));
    let udb = udb.to_str().unwrap();
    let vargs = ["--udbinfo", udb, "--quiet", "--threads", "1"];

    // both report the same header
    let ours = Command::new(vsearch::tool("udbinfo"))
        .args(["-i", udb])
        .output()
        .unwrap();
    let log = tmp("udbinfo-bench.log");
    let _ = std::fs::remove_file(&log);
    let mut args = vargs.to_vec();
    args.extend(["--log", log.to_str().unwrap()]);
    vsearch::run(&vsearch, &args);
    let log = std::fs::read_to_string(&log).unwrap();
    let theirs: String = log
        .lines()
        .skip_while(|l| !l.contains("Seqs"))
        .take(8)
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(String::from_utf8(ours.stdout).unwrap(), theirs);

    let mut group = c.benchmark_group("udbinfo/binary");
    group.sample_size(20);
    vsearch::bench(&mut group, "vsearch", &vsearch, &vargs);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
