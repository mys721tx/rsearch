//! Times `vsearch --fastq_chars` and its rsearch pipeline on the file of
//! rsearch's fqchar benchmark, with the default tail length and a long one
//! (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/fqchar/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let input = input_file();
    let fq = input.to_str().unwrap();

    let mut group = c.benchmark_group("fqchar/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for tail in ["4", "15"] {
        // the report goes to stderr, discarded with the progress lines
        let args = ["--fastq_chars", fq, "--fastq_tail", tail];
        vsearch::bench(&mut group, format!("vsearch/tail_{tail}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
