//! Benchmarks of orienting sequences.
//!
//! `binary` runs the dbornt program end to end on 50,000 reads of 250 bases
//! with one error in fifty. The reads are cut from 50 templates of 1,400
//! bases, half of them reverse-complemented. The database holds 5,000
//! mutated copies of the templates, read from a FASTA file and from mkudb's
//! UDB file of it. The word length is the default of 12. It runs at each
//! thread count of `data::thread_counts` (one and one per core, or
//! `RSEARCH_BENCH_THREADS`). compare/benches/dbornt times vsearch, which
//! orients on one thread, on the same files, after checking that both
//! write the same table.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let (fa, dbfa, udb) = input_files(env!("CARGO_BIN_EXE_mkudb"));
    let fa = fa.to_str().unwrap();
    let dbornt = env!("CARGO_BIN_EXE_dbornt");
    let mut group = c.benchmark_group(data::group("dbornt/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    for (name, db) in [("fasta", &dbfa), ("udb", &udb)] {
        for threads in data::thread_counts() {
            let t = threads.to_string();
            let rargs = [
                "-i",
                fa,
                "--db",
                db.to_str().unwrap(),
                "-o",
                "/dev/null",
                "-t",
                &t,
            ];
            group.bench_function(format!("dbornt/{name}/t{t}"), |b| {
                b.iter(|| run(dbornt, &rargs))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
