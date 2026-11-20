//! Benchmarks of joining read pairs with a pad.
//!
//! `binary` runs the joinpa program end to end on a file of forward FASTQ
//! reads and a file of their mates. Set `VSEARCH` to a vsearch executable
//! to time `vsearch --fastq_join` on the same files for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench joinpa
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// Pairs of reads of `LENGTH` nucleotides.
const PAIRS: usize = 100_000;
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

/// A FASTQ read with an Illumina-like header ending in /1 or /2.
fn read(rng: &mut Rng, i: usize, mate: u8) -> String {
    let seq: String = (0..LENGTH)
        .map(|_| char::from(b"ACGTACGTACGTN"[rng.below(13)]))
        .collect();
    let qual: String = (0..LENGTH)
        .map(|_| char::from(35 + rng.below(39) as u8))
        .collect();
    format!("@M00123:45:000000000-ABCDE:1:1101:{i}:1000/{mate}\n{seq}\n+\n{qual}\n")
}

fn input_files() -> (PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let (mut forward, mut reverse) = (String::new(), String::new());
    for i in 0..PAIRS {
        forward.push_str(&read(&mut rng, i, 1));
        reverse.push_str(&read(&mut rng, i, 2));
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let paths = (dir.join("joinpa-bench-1.fq"), dir.join("joinpa-bench-2.fq"));
    std::fs::write(&paths.0, forward).unwrap();
    std::fs::write(&paths.1, reverse).unwrap();
    paths
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
    let (forward, reverse) = input_files();
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());
    let joinpa = env!("CARGO_BIN_EXE_joinpa");

    let mut group = c.benchmark_group("joinpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    let args = ["-i", fwd, "--reverse", rev, "-o", "/dev/null"];
    group.bench_function("joinpa", |b| b.iter(|| run(joinpa, &args)));
    if let Ok(vsearch) = std::env::var("VSEARCH") {
        let args = [
            "--fastq_join",
            fwd,
            "--reverse",
            rev,
            "--quiet",
            "--fastqout",
            "/dev/null",
        ];
        translated::bench(&mut group, "vsearch", &args);
        group.bench_function("vsearch", |b| b.iter(|| run(&vsearch, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
