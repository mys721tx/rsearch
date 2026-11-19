//! Benchmarks of shuffling.
//!
//! `shuffle` measures `rsearch::random::shuffle` on a deck of a million;
//! `binary` runs the shufl program end to end on files of FASTA and FASTQ
//! reads. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --shuffle` on the same files for comparison (vsearch writes
//! FASTQ input as FASTA):
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench shufl
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::random::{Mt19937_64, shuffle};

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

/// Generates reads as FASTQ and as FASTA.
fn generate() -> (String, String) {
    let mut rng = Rng(0x5eed);
    let (mut fastq, mut fasta) = (String::new(), String::new());
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i}\n{seq}\n"));
    }
    (fastq, fasta)
}

fn bench_shuffle(c: &mut Criterion) {
    let mut group = c.benchmark_group("shufl/shuffle");
    const DECK: usize = 1_000_000;
    group.throughput(Throughput::Elements(DECK as u64));
    let mut deck: Vec<u32> = (0..DECK as u32).collect();
    group.bench_function("mt19937_64", |b| {
        b.iter(|| {
            shuffle(&mut deck, &mut Mt19937_64::new(42));
            black_box(&deck);
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
        fastq: dir.join("shufl-bench.fq"),
        fasta: dir.join("shufl-bench.fa"),
    };
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.fasta, fasta).unwrap();
    inputs
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let shufl = env!("CARGO_BIN_EXE_shufl");
    let vsearch = std::env::var("VSEARCH").ok();
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("shufl/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, file) in [("fasta", fa), ("fastq", fq)] {
        group.bench_function(format!("shufl/{name}"), |b| {
            b.iter(|| run(shufl, &["--seed", "42", "-i", file, "-o", "/dev/null"]))
        });
        if let Some(vsearch) = &vsearch {
            let args = [
                "--shuffle",
                file,
                "--quiet",
                "--notrunclabels",
                "--minseqlength",
                "0",
                "--fasta_width",
                "0",
                "--randseed",
                "42",
                "--output",
                "/dev/null",
            ];
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_shuffle, bench_binary);
criterion_main!(benches);
