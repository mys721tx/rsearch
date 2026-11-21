//! Benchmarks of dereplication.
//!
//! `cluster` measures the in-memory work (`Clusters::add_fastx` and
//! sorting) on parsed records; `binary` runs the derep program end to end
//! on files. The cases are FASTA and FASTQ, with and without `--strand`,
//! without sizes, in input order, Phred+64, the best quality, UC and
//! tab-separated tables, identical labels and prefixes.
//! compare/benches/derep times vsearch on the same files.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::cluster::{Cluster, Clusters, sort_by_size};
use rsearch::fastx::{Reader, Record};

mod input;
use input::*;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn dereplicate(reads: &[Record], strand: bool) -> Vec<Cluster> {
    let mut clusters = Clusters::default();
    for read in reads {
        clusters.add_fastx(read, strand);
    }
    let mut clusters = clusters.into_clusters();
    sort_by_size(&mut clusters);
    clusters
}

fn bench_cluster(c: &mut Criterion) {
    if skip() {
        return;
    }
    let fasta_reads = parse(&text(false));
    let fastq_reads = parse(&text(true));

    let mut group = c.benchmark_group(data::group("derep/cluster"));
    group.throughput(Throughput::Elements(fasta_reads.len() as u64));
    for (name, reads, strand) in [
        ("fasta", &fasta_reads, false),
        ("fasta_strand", &fasta_reads, true),
        ("fastq", &fastq_reads, false),
        ("fastq_strand", &fastq_reads, true),
    ] {
        group.bench_function(name, |b| b.iter(|| black_box(dereplicate(reads, strand))));
    }
    group.finish();
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

/// A benchmark case: derep's arguments.
struct Case<'a> {
    name: &'static str,
    derep: Vec<&'a str>,
}

fn cases<'a>(inputs: &'a Inputs) -> Vec<Case<'a>> {
    let fa = inputs.fasta.to_str().unwrap();
    let fq = inputs.fastq.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let case = |name, file, derep: &[&'a str]| Case {
        name,
        derep: [&["-i", file, "-o", "/dev/null"], derep].concat(),
    };
    vec![
        case("fasta", fa, &[]),
        case("fasta_strand", fa, &["--strand"]),
        case("fastq", fq, &[]),
        case("fastq_strand", fq, &["--strand"]),
        case("unsized", fa, &["--no-sizein", "--no-sizeout"]),
        case("unsorted", fa, &["--unsorted"]),
        case("phred64", fq64, &["--ascii", "64", "--ascii-out", "64"]),
        case("best_qual", fq, &["--best-qual"]),
        case(
            "uc",
            fa,
            &["--strand", "--uc", "/dev/null", "--tabbedout", "/dev/null"],
        ),
        case("same_label", fa, &["--same-label"]),
        case("prefix", fa, &["--prefix"]),
    ]
}

fn bench_binary(c: &mut Criterion) {
    if skip() {
        return;
    }
    let inputs = input_files();
    let derep = env!("CARGO_BIN_EXE_derep");

    let mut group = c.benchmark_group(data::group("derep/binary"));
    group.throughput(Throughput::Elements(data::count(&inputs.fastq, READS)));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("derep/{}", case.name), |b| {
            b.iter(|| run(derep, &case.derep))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_cluster, bench_binary);
criterion_main!(benches);
