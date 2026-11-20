//! Benchmarks of sorting.
//!
//! `binary` runs the sortby program end to end on files of FASTA and FASTQ
//! reads of varied lengths and abundances, by length and by abundance. Set
//! `VSEARCH` to a vsearch executable to time `vsearch --sortbylength` and
//! `--sortbysize` on the same files for comparison (vsearch writes FASTQ
//! input as FASTA):
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench sortby
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// Reads of 150 to 300 nucleotides.
const READS: usize = 100_000;

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

/// Generates reads as FASTQ and as FASTA, with lengths from 150 to 300 and
/// abundances from 1 to 1000, so that both keys have ties.
fn generate() -> (String, String) {
    let mut rng = Rng(0x5eed);
    let (mut fastq, mut fasta) = (String::new(), String::new());
    for i in 0..READS {
        let length = 150 + rng.below(151);
        let seq: String = (0..length)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..length)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        let size = 1 + rng.below(1000);
        fastq.push_str(&format!("@read{i};size={size}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i};size={size}\n{seq}\n"));
    }
    (fastq, fasta)
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
        fastq: dir.join("sortby-bench.fq"),
        fasta: dir.join("sortby-bench.fa"),
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
    let sortby = env!("CARGO_BIN_EXE_sortby");
    let vsearch = std::env::var("VSEARCH").ok();
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("sortby/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (format, file) in [("fasta", fa), ("fastq", fq)] {
        for key in ["length", "size"] {
            let name = format!("{key}_{format}");
            let flag = format!("--{key}");
            group.bench_function(format!("sortby/{name}"), |b| {
                b.iter(|| run(sortby, &[&flag, "-i", file, "-o", "/dev/null"]))
            });
            if let Some(vsearch) = &vsearch {
                let command = format!("--sortby{key}");
                let args = [
                    command.as_str(),
                    file,
                    "--quiet",
                    "--notrunclabels",
                    "--minseqlength",
                    "0",
                    "--fasta_width",
                    "0",
                    "--output",
                    "/dev/null",
                ];
                translated::bench(&mut group, format!("vsearch/{name}"), &args);
                group.bench_function(format!("vsearch/{name}"), |b| {
                    b.iter(|| run(vsearch, &args))
                });
            }
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
