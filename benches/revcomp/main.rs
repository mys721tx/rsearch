//! Benchmarks of reverse complementing.
//!
//! `complement` measures `reverse_complement_keeping_case` on parsed
//! reads; `binary` runs the revcomp program end to end on files of FASTQ
//! and FASTA reads. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --fastx_revcomp` on the same files for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench revcomp
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Reader, Record};

/// Reads of `LENGTH` nucleotides.
const READS: usize = 100_000;
const LENGTH: usize = 250;

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

/// Generates reads as FASTQ and as FASTA, mostly ACGT with some lower case
/// and ambiguity codes.
fn generate() -> (String, String) {
    const CODES: &[u8] = b"ACGTACGTACGTACGTACGTACGTacgtNRY";
    let mut rng = Rng(0x5eed);
    let (mut fastq, mut fasta) = (String::new(), String::new());
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(CODES[rng.below(CODES.len())]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i}\n{seq}\n"));
    }
    (fastq, fasta)
}

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_complement(c: &mut Criterion) {
    let reads = parse(&generate().0);
    let mut group = c.benchmark_group("revcomp/complement");
    group.throughput(Throughput::Elements(READS as u64));
    let mut dest = Vec::new();
    group.bench_function("keeping_case", |b| {
        b.iter(|| {
            for r in &reads {
                reverse_complement_keeping_case(&mut dest, r.seq());
                black_box(&dest);
            }
        })
    });
    group.finish();
}

/// The generated input files, written once into the target directory.
struct Inputs {
    fastq: PathBuf,
    fasta: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fastq, fasta) = generate();
    let inputs = Inputs {
        fastq: dir.join("revcomp-bench.fq"),
        fasta: dir.join("revcomp-bench.fa"),
    };
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.fasta, fasta).unwrap();
    inputs
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let revcomp = env!("CARGO_BIN_EXE_revcomp");
    let vsearch = std::env::var("VSEARCH").ok();
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("revcomp/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, file, output) in [
        ("fastq", fq, &["--fastqout", "/dev/null"][..]),
        (
            "fasta",
            fa,
            &["--fastaout", "/dev/null", "--fasta_width", "0"],
        ),
    ] {
        group.bench_function(format!("revcomp/{name}"), |b| {
            b.iter(|| run(revcomp, &["-i", file, "-o", "/dev/null"]))
        });
        if let Some(vsearch) = &vsearch {
            let args = [&["--fastx_revcomp", file, "--quiet"][..], output].concat();
            translated::bench(&mut group, format!("vsearch/{name}"), &args);
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_complement, bench_binary);
criterion_main!(benches);
