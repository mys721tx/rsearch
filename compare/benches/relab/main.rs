//! Times `vsearch --fastx_filter` with relab's label options and its
//! rsearch pipelines on the files of rsearch's relab benchmark (see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/relab/input.rs"]
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

    let mut group = c.benchmark_group(data::group("relab/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for (name, file, options) in [
        ("copy", fq, &[][..]),
        ("prefix", fq, &["--relabel", "R", "--sizeout"]),
        ("self", fq, &["--relabel_self"]),
        ("sha1", fq, &["--relabel_sha1"]),
        ("md5_keep", fq, &["--relabel_md5", "--relabel_keep"]),
        ("strip", fq, &["--xsize", "--xee", "--xlength"]),
        (
            "annotate",
            fq,
            &[
                "--sizeout",
                "--lengthout",
                "--sample",
                "S1",
                "--label_suffix",
                ";",
            ],
        ),
        ("ee", fq, &["--eeout"]),
        ("fasta_prefix", fa, &["--relabel", "R", "--sizeout"]),
        ("fasta_strip", fa, &["--xsize", "--xee", "--xlength"]),
    ] {
        let output: &[&str] = if file == fa {
            &["--fastaout", "/dev/null", "--fasta_width", "0"]
        } else {
            &["--fastqout", "/dev/null"]
        };
        let args = [&["--fastx_filter", file, "--quiet"], output, options].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
