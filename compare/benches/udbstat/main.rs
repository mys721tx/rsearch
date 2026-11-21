//! Times `vsearch --udbstats --log` (vsearch analyses nothing without a
//! log) and its rsearch pipeline on the UDB file of rsearch's udbstat
//! benchmark. It first checks that udbstat reports the statistics vsearch
//! logs (see `common/vsearch.rs`).

use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};

#[path = "../../../benches/udbstat/input.rs"]
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

    // both report the same statistics
    let log = tmp("udbstat-bench.log");
    let vargs = [
        "--udbstats",
        udb,
        "--quiet",
        "--threads",
        "1",
        "--log",
        log.to_str().unwrap(),
    ];
    let ours = Command::new(vsearch::tool("udbstat"))
        .args(["-i", udb])
        .output()
        .unwrap();
    let _ = std::fs::remove_file(&log);
    vsearch::run(&vsearch, &vargs);
    let text = std::fs::read_to_string(&log).unwrap();
    // the report: past the log's three opening lines, up to its blank line
    // and three closing ones
    let lines: Vec<&str> = text.lines().collect();
    let theirs: String = lines[3..lines.len() - 4]
        .iter()
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(String::from_utf8(ours.stdout).unwrap(), theirs);

    let mut group = c.benchmark_group("udbstat/binary");
    group.sample_size(10);
    vsearch::bench(&mut group, "vsearch", &vsearch, &vargs);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
