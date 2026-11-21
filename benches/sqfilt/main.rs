//! Benchmarks of read filtering.
//!
//! `filter` measures `Filter::apply` with filtering settings on parsed
//! reads; `binary` runs the sqfilt program end to end on files. The cases
//! are no filter, length, Ns, expected errors and their rate, lowest score,
//! abundance, all of them together, pairs, the share of upper-case residues
//! and FASTA. compare/benches/sqfilt times vsearch with the same filters on
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

fn bench_filter(c: &mut Criterion) {
    let reads = parse(&generate("", 0x5eed));
    let mut group = c.benchmark_group("sqfilt/filter");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, filter) in [
        (
            "max_ee",
            Filter {
                max_ee: Some(1.0),
                ..Filter::default()
            },
        ),
        (
            "min_qual",
            Filter {
                min_qual: 10,
                ..Filter::default()
            },
        ),
        (
            "all",
            Filter {
                min_len: 100,
                max_ns: Some(2),
                max_ee: Some(2.0),
                max_rate: Some(0.01),
                min_qual: 5,
                ..Filter::default()
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

/// A benchmark case: sqfilt's arguments.
struct Case<'a> {
    name: &'static str,
    sqfilt: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let rev = inputs.reverse.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let null = "/dev/null";
    let case = |name, file, sqfilt: &[&'static str]| Case {
        name,
        sqfilt: [&["-i", file, "-o", null, "--discarded", null], sqfilt].concat(),
    };
    vec![
        case("none", fq, &[]),
        case("min_len", fq, &["--min-len", "240"]),
        case("max_ns", fq, &["--max-ns", "10"]),
        case("max_ee", fq, &["--max-ee", "1"]),
        case("max_ee_rate", fq, &["--max-ee-rate", "0.005"]),
        case("min_qual", fq, &["--min-qual", "10"]),
        case("size", fq, &["--min-size", "2", "--max-size", "4"]),
        case(
            "all",
            fq,
            &[
                "--min-len",
                "100",
                "--max-ns",
                "20",
                "--max-ee",
                "3",
                "--min-qual",
                "5",
                "--min-size",
                "2",
            ],
        ),
        Case {
            name: "pairs",
            sqfilt: vec![
                "-i",
                fq,
                "--reverse",
                rev,
                "-o",
                null,
                "--rev-out",
                null,
                "--discarded",
                null,
                "--rev-discarded",
                null,
                "--max-ee",
                "1",
            ],
        },
        Case {
            name: "upper_pct",
            sqfilt: vec!["-i", fq, "-o", null, "--min-upper-pct", "90"],
        },
        case("fasta", fa, &["--min-len", "100", "--max-ns", "10"]),
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let sqfilt = env!("CARGO_BIN_EXE_sqfilt");

    let mut group = c.benchmark_group("sqfilt/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("sqfilt/{}", case.name), |b| {
            b.iter(|| run(sqfilt, &case.sqfilt))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_filter, bench_binary);
criterion_main!(benches);
