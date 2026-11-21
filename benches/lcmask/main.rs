//! Benchmarks of DUST masking.
//!
//! `dust` measures `rsearch::dust::dust` on parsed reads, soft and hard;
//! `binary` runs the lcmask program end to end on a FASTQ file, soft and
//! hard. The reads are random with a repeat, a homopolymer or nothing
//! inserted, so that some windows are masked and most are searched in
//! full. compare/benches/lcmask times vsearch on the same file.

use std::hint::black_box;
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::dust::{Style, dust};
use rsearch::fastx::{Reader, Record};

mod input;
use input::*;

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_dust(c: &mut Criterion) {
    let reads = parse(&generate());
    let mut group = c.benchmark_group("lcmask/dust");
    group.throughput(Throughput::Elements(READS as u64));
    let mut seq = Vec::new();
    for (name, style) in [("soft", Style::Soft), ("hard", Style::Hard)] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    seq.clear();
                    seq.extend_from_slice(r.seq());
                    dust(&mut seq, style);
                    black_box(&seq);
                }
            })
        });
    }
    group.finish();
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fq = input.to_str().unwrap();
    let lcmask = env!("CARGO_BIN_EXE_lcmask");

    let mut group = c.benchmark_group("lcmask/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, hard) in [("soft", &[][..]), ("hard", &["--hard"])] {
        let args = [&["-i", fq, "-o", "/dev/null"][..], hard].concat();
        group.bench_function(format!("lcmask/{name}"), |b| b.iter(|| run(lcmask, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_dust, bench_binary);
criterion_main!(benches);
