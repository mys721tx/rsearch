//! Times `vsearch --fastq_eestats2` and its rsearch pipeline on the file of
//! rsearch's eegrid benchmark. It runs with the default cutoffs, and with a
//! cutoff at every length and five expected-error cutoffs (see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/eegrid/input.rs"]
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

    let mut group = c.benchmark_group(data::group("eegrid/binary"));
    group.throughput(Throughput::Elements(data::count(&input, READS)));
    group.sample_size(10);
    for (name, cutoffs) in [
        ("default", &[][..]),
        (
            "fine",
            &["--length_cutoffs", "1,*,1", "--ee_cutoffs", "0.1,0.5,1,2,5"],
        ),
    ] {
        let args = [
            &["--fastq_eestats2", fq, "--quiet", "--output", "/dev/null"][..],
            cutoffs,
        ]
        .concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
