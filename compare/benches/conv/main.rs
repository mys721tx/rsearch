//! Times `vsearch --fastq_convert`, `--fasta2fastq`, `--sff_convert` and
//! `--udb2fasta` and their rsearch pipelines on the files of rsearch's conv
//! benchmark. For FASTA output it times `--fastx_subsample --sample_pct
//! 100`, which writes every read. It first checks that conv writes the
//! FASTA `--udb2fasta` writes (see `common/vsearch.rs`).

use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

#[path = "../../../benches/conv/input.rs"]
mod input;
#[path = "../common/vsearch.rs"]
mod vsearch;
use input::*;

/// The vsearch commands, by name.
fn cases(inputs: &Inputs) -> Vec<(&'static str, Vec<&str>)> {
    let fq = inputs.phred33.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let solexa = inputs.solexa.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let sff = inputs.sff.to_str().unwrap();
    let udb = inputs.udb.to_str().unwrap();
    let case = |name, file, command, options: &[&'static str]| {
        let args = [
            &[command, file, "--fastqout", "/dev/null", "--quiet"],
            options,
        ]
        .concat();
        (name, args)
    };
    let fasta = |name, file, options: &[&'static str]| {
        let args = [
            &[
                "--fastx_subsample",
                file,
                "--sample_pct",
                "100",
                "--minseqlength",
                "0",
                "--notrunclabels",
                "--fastaout",
                "/dev/null",
                "--quiet",
            ],
            options,
        ]
        .concat();
        (name, args)
    };
    let convert = "--fastq_convert";
    vec![
        case(
            "to_phred64",
            fq,
            convert,
            &["--fastq_asciiout", "64", "--fastq_qmaxout", "62"],
        ),
        case("from_phred64", fq64, convert, &["--fastq_ascii", "64"]),
        case(
            "solexa",
            solexa,
            convert,
            &["--fastq_ascii", "64", "--fastq_solexa"],
        ),
        case(
            "limited",
            fq,
            convert,
            &["--fastq_qminout", "2", "--fastq_qmaxout", "40"],
        ),
        case("fasta", fa, "--fasta2fastq", &[]),
        case("sff", sff, "--sff_convert", &[]),
        case("sff_clip", sff, "--sff_convert", &["--sff_clip"]),
        fasta("fastq_to_fasta", fq, &[]),
        fasta("fastq_to_fasta_unwrapped", fq, &["--fasta_width", "0"]),
        fasta("fasta_rewrap", fa, &[]),
        (
            "udb_to_fasta",
            vec!["--udb2fasta", udb, "--output", "/dev/null", "--quiet"],
        ),
    ]
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let Some(vsearch) = vsearch::vsearch() else {
        return;
    };
    let inputs = input_files(&vsearch::tool("mkudb"));

    // conv writes the FASTA vsearch --udb2fasta writes
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (r, v) = (
        data::output(dir.join("conv-bench-r.fa")),
        data::output(dir.join("conv-bench-v.fa")),
    );
    let udb = inputs.udb.to_str().unwrap();
    vsearch::run(
        &vsearch::tool("conv"),
        &["-i", udb, "--fasta", "-o", r.to_str().unwrap()],
    );
    vsearch::run(
        &vsearch,
        &[
            "--udb2fasta",
            udb,
            "--quiet",
            "--output",
            v.to_str().unwrap(),
        ],
    );
    assert!(
        std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
        "conv and vsearch wrote different FASTA"
    );

    let mut group = c.benchmark_group(data::group("conv/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.phred33, READS)));
    group.sample_size(10);
    for (name, args) in cases(&inputs) {
        vsearch::bench(&mut group, format!("vsearch/{name}"), &vsearch, &args);
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
