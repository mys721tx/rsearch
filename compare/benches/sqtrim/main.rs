//! Times `vsearch --fastx_filter` with sqtrim's trimming options and its
//! rsearch pipelines on the files of rsearch's sqtrim benchmark. It writes
//! the discarded reads too, since sqtrim writes every read (see
//! `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/sqtrim/input.rs"]
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
    let fq64 = inputs.phred64.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("sqtrim/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, file, options) in [
        (
            "strip",
            fq,
            &["--fastq_stripleft", "10", "--fastq_stripright", "10"][..],
        ),
        ("len", fq, &["--fastq_trunclen_keep", "200"]),
        ("qual", fq, &["--fastq_truncqual", "10"]),
        ("max_ee", fq, &["--fastq_truncee", "1"]),
        ("max_ee_rate", fq, &["--fastq_truncee_rate", "0.01"]),
        (
            "all",
            fq,
            &[
                "--fastq_stripleft",
                "5",
                "--fastq_trunclen_keep",
                "200",
                "--fastq_truncqual",
                "5",
                "--fastq_truncee",
                "2",
            ],
        ),
        (
            "phred64",
            fq64,
            &["--fastq_ascii", "64", "--fastq_truncqual", "10"],
        ),
        (
            "fasta",
            fa,
            &["--fastq_stripleft", "10", "--fastq_trunclen_keep", "200"],
        ),
    ] {
        let output: &[&str] = if file == fa {
            &[
                "--fastaout",
                "/dev/null",
                "--fastaout_discarded",
                "/dev/null",
                "--fasta_width",
                "0",
            ]
        } else {
            &[
                "--fastqout",
                "/dev/null",
                "--fastqout_discarded",
                "/dev/null",
            ]
        };
        let args = [&["--fastx_filter", file, "--quiet"], output, options].concat();
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
