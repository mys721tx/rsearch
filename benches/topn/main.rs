//! Benchmarks of keeping the first records.
//!
//! `binary` runs the topn program end to end on files of FASTQ and FASTA
//! reads, keeping 10 records or all of them. Either way it reads and checks
//! the whole input. vsearch has no command of its own for this, so with
//! `VSEARCH` set to a vsearch executable, the pipeline `derep | topn 10` is
//! timed against `vsearch --fastx_uniques --topn 10` on the same FASTA file:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench topn
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

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

/// Generates reads as FASTQ and as FASTA, drawn from 5000 distinct
/// sequences so that dereplication has duplicates to merge.
fn generate() -> (String, String) {
    let mut rng = Rng(0x5eed);
    let distinct: Vec<String> = (0..5000)
        .map(|_| {
            (0..LENGTH)
                .map(|_| char::from(b"ACGT"[rng.below(4)]))
                .collect()
        })
        .collect();
    let (mut fastq, mut fasta) = (String::new(), String::new());
    for i in 0..READS {
        let seq = &distinct[rng.below(distinct.len())];
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i}\n{seq}\n"));
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
        fastq: dir.join("topn-bench.fq"),
        fasta: dir.join("topn-bench.fa"),
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
    let topn = env!("CARGO_BIN_EXE_topn");
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("topn/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, file, count) in [
        ("fastq_10", fq, "10"),
        ("fastq_all", fq, "1000000"),
        ("fasta_10", fa, "10"),
        ("fasta_all", fa, "1000000"),
    ] {
        group.bench_function(format!("topn/{name}"), |b| {
            b.iter(|| run(topn, &["-i", file, "-o", "/dev/null", count]))
        });
    }

    if let Ok(vsearch) = std::env::var("VSEARCH") {
        let derep = env!("CARGO_BIN_EXE_derep");
        let pipeline = format!("{derep} -i {fa} | {topn} 10 > /dev/null");
        group.bench_function("pipeline/derep_topn_10", |b| {
            b.iter(|| run("sh", &["-c", &pipeline]))
        });
        group.bench_function("vsearch/fastx_uniques_topn_10", |b| {
            b.iter(|| {
                run(
                    &vsearch,
                    &[
                        "--fastx_uniques",
                        fa,
                        "--quiet",
                        "--sizein",
                        "--sizeout",
                        "--notrunclabels",
                        "--minseqlength",
                        "0",
                        "--fasta_width",
                        "0",
                        "--topn",
                        "10",
                        "--fastaout",
                        "/dev/null",
                    ],
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
