//! Benchmarks of the grid of reads kept by length and expected errors.
//!
//! `binary` runs the eegrid program end to end on a file of FASTQ reads of
//! varied lengths. It runs with the default cutoffs, and with a cutoff at
//! every length and five expected-error cutoffs. compare/benches/eegrid
//! times vsearch on the same file.

use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

mod input;
use input::*;

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fq = input.to_str().unwrap();
    let eegrid = env!("CARGO_BIN_EXE_eegrid");

    let mut group = c.benchmark_group("eegrid/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, cutoffs) in [
        ("default", &[][..]),
        (
            "fine",
            &[
                "--min-len",
                "1",
                "--len-step",
                "1",
                "--max-ee",
                "0.1,0.5,1,2,5",
            ],
        ),
    ] {
        let args = [&["-i", fq, "-o", "/dev/null"][..], cutoffs].concat();
        group.bench_function(format!("eegrid/{name}"), |b| b.iter(|| run(eegrid, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
