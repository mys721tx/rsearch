//! Times `vsearch --shuffle` and its rsearch pipelines on the files of
//! rsearch's shufl benchmark (vsearch writes FASTQ input as FASTA; see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/shufl/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let inputs = input_files();
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("shufl/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, file) in [("fasta", fa), ("fastq", fq)] {
        let args = [
            "--shuffle",
            file,
            "--quiet",
            "--notrunclabels",
            "--minseqlength",
            "0",
            "--fasta_width",
            "0",
            "--randseed",
            "42",
            "--output",
            "/dev/null",
        ];
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
