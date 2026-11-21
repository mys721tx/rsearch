//! Benchmarks of label rewriting.
//!
//! `format` measures `Label::format` on parsed reads; `binary` runs the
//! relab program end to end on FASTQ and FASTA files. The cases are
//! copying, relabelling by prefix, sequence, SHA1 and MD5, removing and
//! adding annotations, and expected errors. compare/benches/relab times
//! vsearch with the same label options on the same files.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::label::{Label, Relabel, abundance};

mod input;
use input::*;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_format(c: &mut Criterion) {
    if skip() {
        return;
    }
    let reads = parse(&text());
    let mut group = c.benchmark_group(data::group("relab/format"));
    group.throughput(Throughput::Elements(reads.len() as u64));
    for (name, label) in [
        (
            "prefix",
            Label {
                relabel: Relabel::Prefix("R".into()),
                sizeout: true,
                ..Label::default()
            },
        ),
        (
            "sha1",
            Label {
                relabel: Relabel::Sha1,
                ..Label::default()
            },
        ),
        (
            "strip",
            Label {
                xsize: true,
                xee: true,
                xlength: true,
                ..Label::default()
            },
        ),
        (
            "annotate",
            Label {
                sizeout: true,
                lengthout: true,
                sample: Some("S1".into()),
                ..Label::default()
            },
        ),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for (i, r) in reads.iter().enumerate() {
                    let size = abundance(r.header());
                    black_box(label.format(r.header(), r.seq(), size, i as u64 + 1, None));
                }
            })
        });
    }
    group.finish();
}

/// A benchmark case: relab's arguments.
struct Case<'a> {
    name: &'static str,
    relab: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let case = |name, file, relab: &[&'static str]| Case {
        name,
        relab: [&["-i", file, "-o", "/dev/null"], relab].concat(),
    };
    vec![
        case("copy", fq, &[]),
        case("prefix", fq, &["--relabel", "R", "--sizeout"]),
        case("self", fq, &["--relabel-self"]),
        case("sha1", fq, &["--relabel-sha1"]),
        case("md5_keep", fq, &["--relabel-md5", "--relabel-keep"]),
        case("strip", fq, &["--xsize", "--xee", "--xlength"]),
        case(
            "annotate",
            fq,
            &[
                "--sizeout",
                "--lengthout",
                "--sample",
                "S1",
                "--label-suffix",
                ";",
            ],
        ),
        case("ee", fq, &["--ee"]),
        case("fasta_prefix", fa, &["--relabel", "R", "--sizeout"]),
        case("fasta_strip", fa, &["--xsize", "--xee", "--xlength"]),
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let inputs = input_files();
    let relab = env!("CARGO_BIN_EXE_relab");

    let mut group = c.benchmark_group(data::group("relab/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("relab/{}", case.name), |b| {
            b.iter(|| run(relab, &case.relab))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_format, bench_binary);
criterion_main!(benches);
