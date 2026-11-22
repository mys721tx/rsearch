//! Benchmarks of quality conversion.
//!
//! `converter` measures `Converter::convert` on parsed reads; `binary` runs
//! the conv program end to end on files. The cases are Phred+33 to Phred+64
//! and back, Solexa scores, limited output scores, FASTA to FASTQ, SFF to
//! FASTQ with and without clipping, FASTQ and FASTA to wrapped or unwrapped
//! FASTA, and UDB (written by mkudb from the FASTA reads) to FASTA.
//! compare/benches/conv times vsearch on the same files.
//!
//! `simd` times conv from FASTA and FASTQ to FASTA at each SIMD level this
//! processor has, chosen with `RSEARCH_SIMD` or listed in `RSEARCH_BENCH_SIMD`,
//! for the reader. It sets the scalar code the compiler vectorises against the
//! hand-written vector code.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::quality::Conversion;

mod input;
use input::*;
#[path = "../common/levels.rs"]
mod levels;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_converter(c: &mut Criterion) {
    if skip() {
        return;
    }
    let reads = parse(&text());
    let mut group = c.benchmark_group(data::group("conv/converter"));
    group.throughput(Throughput::Elements(reads.len() as u64));
    for (name, conversion) in [
        (
            "phred64",
            Conversion {
                offset_out: 64,
                max_out: 62,
                ..Conversion::default()
            },
        ),
        (
            "limited",
            Conversion {
                min_out: 2,
                max_out: 40,
                ..Conversion::default()
            },
        ),
    ] {
        let converter = conversion.converter();
        let mut out = Vec::new();
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    converter.convert(r.qual(), &mut out).unwrap();
                    black_box(&out);
                }
            })
        });
    }
    group.finish();
}

/// A benchmark case: conv's arguments.
struct Case<'a> {
    name: &'static str,
    conv: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.phred33.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let solexa = inputs.solexa.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let sff = inputs.sff.to_str().unwrap();
    let udb = inputs.udb.to_str().unwrap();
    let case = |name, file, conv: &[&'static str]| Case {
        name,
        conv: [&["-i", file, "-o", "/dev/null"], conv].concat(),
    };
    let fasta = |name, file, conv: &[&'static str]| Case {
        name,
        conv: [&["-i", file, "-o", "/dev/null", "--fasta"], conv].concat(),
    };
    vec![
        case("to_phred64", fq, &["--ascii-out", "64", "--qmax-out", "62"]),
        case("from_phred64", fq64, &["--ascii", "64"]),
        case("solexa", solexa, &["--ascii", "64", "--solexa"]),
        case("limited", fq, &["--qmin-out", "2", "--qmax-out", "40"]),
        case("fasta", fa, &[]),
        case("sff", sff, &[]),
        case("sff_clip", sff, &["--clip"]),
        fasta("fastq_to_fasta", fq, &[]),
        fasta("fastq_to_fasta_unwrapped", fq, &["--fasta-width", "0"]),
        fasta("fasta_rewrap", fa, &[]),
        fasta("udb_to_fasta", udb, &[]),
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
    let inputs = input_files(env!("CARGO_BIN_EXE_mkudb"));
    let conv = env!("CARGO_BIN_EXE_conv");

    let mut group = c.benchmark_group(data::group("conv/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.phred33, READS)));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("conv/{}", case.name), |b| {
            b.iter(|| run(conv, &case.conv))
        });
    }
    group.finish();
}

/// conv from FASTA and FASTQ to FASTA at each SIMD level this processor has. It
/// mostly reads, and checks each sequence line.
fn bench_simd(c: &mut Criterion) {
    if skip() {
        return;
    }
    let inputs = input_files(env!("CARGO_BIN_EXE_mkudb"));
    let conv = env!("CARGO_BIN_EXE_conv");
    let mut group = c.benchmark_group(data::group("conv/simd"));
    group.throughput(Throughput::Elements(data::count(&inputs.phred33, READS)));
    group.sample_size(10);
    for (name, input) in [("fasta", &inputs.fasta), ("fastq", &inputs.phred33)] {
        let args = ["-i", input.to_str().unwrap(), "--fasta", "-o", "/dev/null"];
        for level in levels::levels() {
            group.bench_function(format!("conv/{name}/{}", level.name()), |b| {
                b.iter(|| levels::run_at(level, conv, &args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_converter, bench_binary, bench_simd);
criterion_main!(benches);
