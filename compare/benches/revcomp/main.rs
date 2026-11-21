//! Times `vsearch --fastx_revcomp` and its rsearch pipelines on the files
//! of rsearch's revcomp benchmark (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/revcomp/input.rs"]
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

    let mut group = c.benchmark_group(data::group("revcomp/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for (name, file, output) in [
        ("fastq", fq, &["--fastqout", "/dev/null"][..]),
        (
            "fasta",
            fa,
            &["--fastaout", "/dev/null", "--fasta_width", "0"],
        ),
    ] {
        let args = [&["--fastx_revcomp", file, "--quiet"][..], output].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
