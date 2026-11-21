//! Times `vsearch --fastx_uniques`, `--derep_smallmem`, `--derep_id` and
//! `--derep_prefix` and their rsearch pipelines on the files of rsearch's
//! derep benchmark. vsearch runs with `--minseqlength 0 --notrunclabels`,
//! as derep keeps every length and whole headers (see `common/vsearch.rs`).

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/derep/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

/// The vsearch commands, by name.
fn cases(inputs: &Inputs) -> Vec<(&'static str, Vec<&str>)> {
    let fa = inputs.fasta.to_str().unwrap();
    let fq = inputs.fastq.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let common = ["--minseqlength", "0", "--notrunclabels", "--quiet"];
    // vsearch `command` on `file`, written to `out`, with `options`
    let case = |name, file, command, out, options: &[&'static str]| {
        let args = [&[command, file, out, "/dev/null"], &common[..], options].concat();
        (name, args)
    };
    let uniques = |name, file, out, options: &[&'static str]| {
        case(name, file, "--fastx_uniques", out, options)
    };
    let (fasta, fastq) = ("--fastaout", "--fastqout");
    let sizes = ["--sizein", "--sizeout"];
    let strand = ["--sizein", "--sizeout", "--strand", "both"];
    vec![
        uniques("fasta", fa, fasta, &sizes),
        uniques("fasta_strand", fa, fasta, &strand),
        uniques("fastq", fq, fastq, &sizes),
        uniques("fastq_strand", fq, fastq, &strand),
        uniques("unsized", fa, fasta, &[]),
        case("unsorted", fa, "--derep_smallmem", fasta, &sizes),
        uniques(
            "phred64",
            fq64,
            fastq,
            &[
                "--sizein",
                "--sizeout",
                "--fastq_ascii",
                "64",
                "--fastq_asciiout",
                "64",
            ],
        ),
        uniques(
            "best_qual",
            fq,
            fastq,
            &["--sizein", "--sizeout", "--fastq_qout_max"],
        ),
        uniques(
            "uc",
            fa,
            fasta,
            &[
                "--sizein",
                "--sizeout",
                "--strand",
                "both",
                "--uc",
                "/dev/null",
                "--tabbedout",
                "/dev/null",
            ],
        ),
        case("same_label", fa, "--derep_id", "--output", &sizes),
        case("prefix", fa, "--derep_prefix", "--output", &sizes),
    ]
}

fn bench_binary(c: &mut Criterion) {
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let inputs = input_files();

    let mut group = c.benchmark_group("derep/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, args) in cases(&inputs) {
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
