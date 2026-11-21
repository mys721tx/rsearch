//! Benchmarks of formatting hits.
//!
//! `binary` runs the hitfmt program end to end on the hits exsrch finds on
//! both strands. The 100,000 reads of 150 bases are drawn from 20,000
//! distinct sequences, a quarter of them reverse-complemented and a tenth
//! with one base changed. The database holds those 20,000 sequences and
//! 5,000 duplicates. It writes user fields, BLAST tabular and alignments.
//! vsearch has no command that only formats hits; compare/benches/hitfmt
//! times the closest one, `vsearch --search_exact` writing the same output,
//! on the same files.

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
    if skip() {
        return;
    }
    let (fa, db) = input_files();
    let (fa, db) = (fa.to_str().unwrap(), db.to_str().unwrap());
    let hits = data::output(tmp("hitfmt-bench.uc"));
    let hits = hits.to_str().unwrap();
    run(
        env!("CARGO_BIN_EXE_exsrch"),
        &["-i", fa, "--db", db, "-o", hits, "--strand", "both"],
    );
    let hitfmt = env!("CARGO_BIN_EXE_hitfmt");
    let fields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+raw+caln+qrow";
    let mut group = c.benchmark_group(data::group("hitfmt/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    for (name, rargs) in [
        ("userout", vec!["--userout", fields]),
        ("blast6", vec!["--blast6"]),
        ("alignments", vec!["--alignments"]),
    ] {
        // exact hits, scored as vsearch's exact search scores them
        let rargs = [
            &[hits, "-i", fa, "--db", db, "-o", "/dev/null", "--exact"][..],
            &rargs,
        ]
        .concat();
        group.bench_function(format!("hitfmt/{name}"), |b| b.iter(|| run(hitfmt, &rargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
