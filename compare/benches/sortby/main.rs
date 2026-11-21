//! Times `vsearch --sortbylength` and `--sortbysize` and their rsearch
//! pipelines on the files of rsearch's sortby benchmark (vsearch writes
//! FASTQ input as FASTA; see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/sortby/input.rs"]
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
    let inputs = input_files();
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group(data::group("sortby/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for (format, file) in [("fasta", fa), ("fastq", fq)] {
        for key in ["length", "size"] {
            let command = format!("--sortby{key}");
            let args = [
                command.as_str(),
                file,
                "--quiet",
                "--notrunclabels",
                "--minseqlength",
                "0",
                "--fasta_width",
                "0",
                "--output",
                "/dev/null",
            ];
            vsearch::bench(
                &mut group,
                format!("vsearch/{key}_{format}"),
                &vsearch,
                &args,
            );
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
