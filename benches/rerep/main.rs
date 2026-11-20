//! Benchmarks of rereplication.
//!
//! `binary` runs the rerep program end to end on files of dereplicated
//! reads, with abundances from 1 to 20: FASTA with and without `;size=1` on
//! the copies, and FASTQ. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --rereplicate` on the FASTA file for comparison (vsearch reads
//! only FASTA):
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench rerep
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// Dereplicated reads of `LENGTH` nucleotides, with abundances from 1 to
/// `MAX_SIZE`.
const READS: usize = 20_000;
const LENGTH: usize = 250;
const MAX_SIZE: usize = 20;

/// A small deterministic generator (xorshift64*). Every run benchmarks the
/// same data without extra dependencies.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Generates dereplicated reads as FASTQ and as FASTA, and returns them
/// with the number of copies they make.
fn generate() -> (String, String, usize) {
    let mut rng = Rng(0x5eed);
    let (mut fastq, mut fasta, mut copies) = (String::new(), String::new(), 0);
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        let size = 1 + rng.below(MAX_SIZE);
        copies += size;
        fastq.push_str(&format!("@read{i};size={size}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i};size={size}\n{seq}\n"));
    }
    (fastq, fasta, copies)
}

/// The generated input files, written once into the target directory, and
/// the number of copies they make.
struct Inputs {
    fastq: PathBuf,
    fasta: PathBuf,
    copies: usize,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fastq, fasta, copies) = generate();
    let inputs = Inputs {
        fastq: dir.join("rerep-bench.fq"),
        fasta: dir.join("rerep-bench.fa"),
        copies,
    };
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.fasta, fasta).unwrap();
    inputs
}

/// A benchmark case: rerep's arguments and the equivalent vsearch ones, if
/// vsearch can run it.
struct Case<'a> {
    name: &'static str,
    rerep: Vec<&'a str>,
    vsearch: Option<Vec<&'a str>>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let null = "/dev/null";
    let vsearch = |extra: &[&'static str]| {
        let common = [
            "--rereplicate",
            fa,
            "--quiet",
            "--notrunclabels",
            "--fasta_width",
            "0",
            "--output",
            null,
        ];
        Some([&common[..], extra].concat())
    };
    vec![
        Case {
            name: "fasta",
            rerep: vec!["-i", fa, "-o", null],
            vsearch: vsearch(&["--sizeout"]),
        },
        Case {
            name: "fasta_no_sizeout",
            rerep: vec!["-i", fa, "-o", null, "--no-sizeout"],
            vsearch: vsearch(&[]),
        },
        Case {
            name: "fastq",
            rerep: vec!["-i", fq, "-o", null],
            vsearch: None,
        },
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let rerep = env!("CARGO_BIN_EXE_rerep");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("rerep/binary");
    group.throughput(Throughput::Elements(inputs.copies as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("rerep/{}", case.name), |b| {
            b.iter(|| run(rerep, &case.rerep))
        });
        if let (Some(vsearch), Some(args)) = (&vsearch, &case.vsearch) {
            translated::bench(&mut group, format!("vsearch/{}", case.name), args);
            group.bench_function(format!("vsearch/{}", case.name), |b| {
                b.iter(|| run(vsearch, args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
