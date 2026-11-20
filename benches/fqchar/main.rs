//! Benchmarks of the residue and score symbol report.
//!
//! `binary` runs the fqchar program end to end on a file of FASTQ reads,
//! with the default tail length and a long one. Set `VSEARCH` to a vsearch
//! executable to time `vsearch --fastq_chars` on the same file for
//! comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench fqchar
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
/// along the read, ending in a run of score 2 now and then.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        let tail = if rng.below(10) == 0 { 20 } else { 0 };
        for j in 0..LENGTH {
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
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("fqchar-bench.fq");
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

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fq = input.to_str().unwrap();
    let fqchar = env!("CARGO_BIN_EXE_fqchar");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("fqchar/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for tail in ["4", "15"] {
        group.bench_function(format!("fqchar/tail_{tail}"), |b| {
            b.iter(|| run(fqchar, &["-i", fq, "-o", "/dev/null", "--tail", tail]))
        });
        if let Some(vsearch) = &vsearch {
            // the report goes to stderr, discarded with the progress lines
            translated::bench(
                &mut group,
                format!("vsearch/tail_{tail}"),
                &["--fastq_chars", fq, "--fastq_tail", tail],
            );
            group.bench_function(format!("vsearch/tail_{tail}"), |b| {
                b.iter(|| run(vsearch, &["--fastq_chars", fq, "--fastq_tail", tail]))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
