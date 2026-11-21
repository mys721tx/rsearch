//! Times `vsearch --fastx_mask`, soft and hard, and its rsearch pipelines
//! on the file of rsearch's lcmask benchmark (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/lcmask/input.rs"]
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

    let mut group = c.benchmark_group("lcmask/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, hardmask) in [("soft", &[][..]), ("hard", &["--hardmask"])] {
        let args = [
            &[
                "--fastx_mask",
                fq,
                "--quiet",
                "--notrunclabels",
                "--minseqlength",
                "0",
                "--threads",
                "1",
                "--fastqout",
                "/dev/null",
            ][..],
            hardmask,
        ]
        .concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
