//! Benchmarks of OTU tables.
//!
//! `binary` runs the otutab program end to end on the hits exsrch finds on
//! both strands. The 100,000 reads of 150 bases come from 20 samples and
//! are drawn from 20,000 distinct sequences of 2,000 OTUs, a quarter of
//! them reverse-complemented and a tenth with one base changed. The
//! database holds those 20,000 sequences and 5,000 duplicates. vsearch has
//! no command that only counts hits; compare/benches/otutab times the
//! closest one, `vsearch --search_exact` writing the same table, on the
//! same files.

use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let hits = tmp("otutab-bench.uc");
    let hits = hits.to_str().unwrap();
    run(
        env!("CARGO_BIN_EXE_exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let otutab = env!("CARGO_BIN_EXE_otutab");
    let mut group = c.benchmark_group("otutab/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, rargs) in [
        ("classic", vec![]),
        ("mothur", vec!["--mothur"]),
        ("biom", vec!["--biom"]),
    ] {
        let rargs = [&[hits, "--db", db, "-o", "/dev/null"][..], &rargs].concat();
        group.bench_function(format!("otutab/{name}"), |b| b.iter(|| run(otutab, &rargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
