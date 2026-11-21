//! Times `vsearch --rereplicate` and its rsearch pipelines on the FASTA
//! file of rsearch's rerep benchmark (vsearch reads only FASTA; see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/rerep/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let inputs = input_files();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("rerep/binary");
    group.throughput(Throughput::Elements(inputs.copies as u64));
    group.sample_size(10);
    for (name, extra) in [("fasta", &["--sizeout"][..]), ("fasta_no_sizeout", &[])] {
        let common = [
            "--rereplicate",
            fa,
            "--quiet",
            "--notrunclabels",
            "--fasta_width",
            "0",
            "--output",
            "/dev/null",
        ];
        let args = [&common[..], extra].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
