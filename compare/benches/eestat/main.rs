//! Times `vsearch --fastq_eestats` and its rsearch pipeline on the file of
//! rsearch's eestat benchmark (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/eestat/input.rs"]
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

    let mut group = c.benchmark_group("eestat/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    let args = ["--fastq_eestats", fq, "--quiet", "--output", "/dev/null"];
    vsearch::bench(&mut group, "vsearch", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
