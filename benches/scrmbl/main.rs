//! Benchmarks of scrambling.
//!
//! `scramble` measures `Scrambler::scramble` on parsed reads, with k of 1,
//! 2 and 6; `binary` runs the scrmbl program end to end on a FASTQ file
//! with the same k. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --scramble` on the same file for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench scrmbl
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::random::{SplitMix64, substream};
use rsearch::scramble::Scrambler;

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

/// Generates FASTQ reads.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
    }
    fastq
}

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_scramble(c: &mut Criterion) {
    let reads = parse(&generate());
    let mut group = c.benchmark_group("scrmbl/scramble");
    group.throughput(Throughput::Elements(READS as u64));
    let mut scrambler = Scrambler::default();
    let mut seq = Vec::new();
    for k in [1, 2, 6] {
        group.bench_function(format!("k{k}"), |b| {
            b.iter(|| {
                for (i, r) in reads.iter().enumerate() {
                    seq.clear();
                    seq.extend_from_slice(r.seq());
                    let mut generator = SplitMix64::new(substream(42, i as u64 + 1));
                    scrambler.scramble(&mut seq, k, &mut generator);
                    black_box(&seq);
                }
            })
        });
    }
    group.finish();
}

fn input_file() -> PathBuf {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("scrmbl-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fq = input.to_str().unwrap();
    let scrmbl = env!("CARGO_BIN_EXE_scrmbl");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("scrmbl/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for k in ["1", "2", "6"] {
        group.bench_function(format!("scrmbl/k{k}"), |b| {
            b.iter(|| {
                run(
                    scrmbl,
                    &["-i", fq, "-o", "/dev/null", "--seed", "42", "--kmer", k],
                )
            })
        });
        if let Some(vsearch) = &vsearch {
            let args = [
                "--scramble",
                fq,
                "--quiet",
                "--notrunclabels",
                "--fastqout",
                "/dev/null",
                "--randseed",
                "42",
                "--scramble_kmer",
                k,
            ];
            group.bench_function(format!("vsearch/k{k}"), |b| b.iter(|| run(vsearch, &args)));
        }
    }
    group.finish();
}

criterion_group!(benches, bench_scramble, bench_binary);
criterion_main!(benches);
