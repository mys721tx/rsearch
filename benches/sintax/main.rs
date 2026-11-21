//! Benchmarks of SINTAX classification.
//!
//! `binary` runs the sintax program end to end: 10,000 reads of 250 bases
//! with one error in fifty, cut from 50 templates of 1,400 bases, half of
//! them reverse-complemented, against a database of 5,000 mutated copies of
//! the templates, each annotated with its template's taxonomy, at the
//! default word length of 8, searching both strands, from a FASTA file and
//! from mkudb's UDB file of it. compare/benches/sintax times vsearch on the
//! same files, after checking that both write the same table.

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let (fa, dbfa, udb) = input_files(env!("CARGO_BIN_EXE_mkudb"));
    let fa = fa.to_str().unwrap();
    let sintax = env!("CARGO_BIN_EXE_sintax");
    let r = data::output(tmp("sintax-bench-r.tsv"));
    let r = r.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("sintax/binary"));
    group.throughput(Throughput::Elements(data::count(fa, READS)));
    group.sample_size(10);
    for (name, db) in [("fasta", &dbfa), ("udb", &udb)] {
        let db = db.to_str().unwrap();
        let rargs = [
            "-i", fa, "--db", db, "-o", r, "--seed", "1", "--strand", "both",
        ];
        group.bench_function(format!("sintax/{name}"), |b| b.iter(|| run(sintax, &rargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
