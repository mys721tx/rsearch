//! Benchmarks of read trimming.
//!
//! `trim` measures `Filter::apply` with trimming settings on parsed reads;
//! `binary` runs the sqtrim program end to end on files. The cases are
//! stripping ends, truncating to a length, at a low score, at an
//! expected-error limit or rate, all of them together, Phred+64 and FASTA.
//! compare/benches/sqtrim times vsearch with the same trimming options on
//! the same files.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::filter::Filter;

mod input;
use input::*;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_trim(c: &mut Criterion) {
    if skip() {
        return;
    }
    let reads = parse(&text());
    let mut group = c.benchmark_group(data::group("sqtrim/trim"));
    group.throughput(Throughput::Elements(reads.len() as u64));
    let trimming = Filter {
        min_len: 0,
        ..Filter::default()
    };
    for (name, filter) in [
        (
            "strip_len",
            Filter {
                strip_left: 10,
                trunc_keep: Some(200),
                ..trimming.clone()
            },
        ),
        (
            "qual",
            Filter {
                trunc_qual: Some(10),
                ..trimming.clone()
            },
        ),
        (
            "max_ee",
            Filter {
                trunc_ee: Some(1.0),
                ..trimming.clone()
            },
        ),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    black_box(filter.apply(r.seq(), Some(r.qual()), 1).unwrap());
                }
            })
        });
    }
    group.finish();
}

/// A benchmark case: sqtrim's arguments.
struct Case<'a> {
    name: &'static str,
    sqtrim: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let case = |name, file, sqtrim: &[&'static str]| Case {
        name,
        sqtrim: [&["-i", file, "-o", "/dev/null"], sqtrim].concat(),
    };
    vec![
        case("strip", fq, &["--strip-left", "10", "--strip-right", "10"]),
        case("len", fq, &["--len", "200"]),
        case("qual", fq, &["--qual", "10"]),
        case("max_ee", fq, &["--max-ee", "1"]),
        case("max_ee_rate", fq, &["--max-ee-rate", "0.01"]),
        case(
            "all",
            fq,
            &[
                "--strip-left",
                "5",
                "--len",
                "200",
                "--qual",
                "5",
                "--max-ee",
                "2",
            ],
        ),
        case("phred64", fq64, &["--ascii", "64", "--qual", "10"]),
        case("fasta", fa, &["--strip-left", "10", "--len", "200"]),
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
    let sqtrim = env!("CARGO_BIN_EXE_sqtrim");

    let mut group = c.benchmark_group(data::group("sqtrim/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("sqtrim/{}", case.name), |b| {
            b.iter(|| run(sqtrim, &case.sqtrim))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_trim, bench_binary);
criterion_main!(benches);
