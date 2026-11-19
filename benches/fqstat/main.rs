//! Benchmarks of the read length and quality report.
//!
//! `binary` runs the fqstat program end to end on a file of FASTQ reads of
//! varied lengths. Set `VSEARCH` to a vsearch executable to time `vsearch
//! --fastq_stats` on the same file for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench fqstat
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

/// Generates Illumina-like reads with scores drifting from high to low
/// along the read. A read ends in a run of score 2 now and then. One read
/// in three is trimmed to a length from 150 to 250.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        let tail = if rng.below(10) == 0 { 20 } else { 0 };
        let length = if i % 3 == 0 {
            LENGTH - rng.below(101)
        } else {
            LENGTH
        };
        for j in 0..length {
            q = (q + rng.below(4) as i64 - 2).clamp(2, 41);
            seq.push(char::from(b"ACGTACGTACGTN"[rng.below(13)]));
            let score = if j + tail >= LENGTH { 2 } else { q };
            qual.push(char::from((score + 33) as u8));
        }
        out.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
    }
    out
}

fn input_file() -> PathBuf {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("fqstat-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}

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
    let fqstat = env!("CARGO_BIN_EXE_fqstat");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("fqstat/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    group.bench_function("fqstat", |b| {
        b.iter(|| run(fqstat, &["-i", fq, "-o", "/dev/null"]))
    });
    if let Some(vsearch) = &vsearch {
        group.bench_function("vsearch", |b| {
            b.iter(|| {
                run(
                    vsearch,
                    &["--fastq_stats", fq, "--quiet", "--log", "/dev/null"],
                )
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
