//! Times `vsearch --fastq_stats` and its rsearch pipeline on the file of
//! rsearch's fqstat benchmark (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/fqstat/input.rs"]
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
    let input = input_file();
    let fq = input.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("fqstat/binary"));
    group.throughput(Throughput::Elements(data::count(&input, READS)));
    group.sample_size(10);
    let args = ["--fastq_stats", fq, "--quiet", "--log", "/dev/null"];
    vsearch::bench(&mut group, "vsearch", &vsearch, &args);
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
