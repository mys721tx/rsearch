//! Benchmarks of pairing the reads of two files by label.
//!
//! `binary` runs the syncpa program end to end on a file of forward FASTQ
//! reads and a file of their mates in shuffled order, with a read in ten
//! missing its mate on each side, writing all four outputs. Set `VSEARCH`
//! to a vsearch executable to time `vsearch --fastx_syncpairs` on the same
//! files for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench syncpa
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// Pairs of reads of `LENGTH` nucleotides.
const PAIRS: usize = 100_000;
const LENGTH: usize = 150;

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

/// A FASTQ read with Illumina-like headers ending in /1 or /2.
fn read(rng: &mut Rng, i: usize, mate: u8) -> String {
    let seq: String = (0..LENGTH)
        .map(|_| char::from(b"ACGT"[rng.below(4)]))
        .collect();
    let qual: String = (0..LENGTH)
        .map(|_| char::from(35 + rng.below(39) as u8))
        .collect();
    format!(
        "@M00123:45:000000000-ABCDE:1:1101:{i}:{}/{mate}\n{seq}\n+\n{qual}\n",
        1000 + i % 9000
    )
}

/// Generates the forward reads in order and their mates shuffled, each
/// side missing a tenth of the mates of the other.
fn generate() -> (String, String) {
    let mut rng = Rng(0x5eed);
    let mut forward = String::new();
    let mut reverse = Vec::new();
    for i in 0..PAIRS {
        let side = rng.below(20);
        if side != 0 {
            forward.push_str(&read(&mut rng, i, 1));
        }
        if side != 1 {
            reverse.push(read(&mut rng, i, 2));
        }
    }
    for i in (1..reverse.len()).rev() {
        let j = rng.below(i + 1);
        reverse.swap(i, j);
    }
    (forward, reverse.concat())
}

fn input_files() -> (PathBuf, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (forward, reverse) = generate();
    let paths = (dir.join("syncpa-bench-1.fq"), dir.join("syncpa-bench-2.fq"));
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

fn bench_binary(c: &mut Criterion) {
    let (forward, reverse) = input_files();
    let (fwd, rev) = (forward.to_str().unwrap(), reverse.to_str().unwrap());
    let syncpa = env!("CARGO_BIN_EXE_syncpa");
    let null = "/dev/null";

    let mut group = c.benchmark_group("syncpa/binary");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.sample_size(10);
    let args = [
        "-i",
        fwd,
        "--reverse",
        rev,
        "-o",
        null,
        "--rev-out",
        null,
        "--orphans",
        null,
        "--rev-orphans",
        null,
    ];
    group.bench_function("syncpa", |b| b.iter(|| run(syncpa, &args)));
    if let Ok(vsearch) = std::env::var("VSEARCH") {
        let args = [
            "--fastx_syncpairs",
            fwd,
            "--reverse",
            rev,
            "--quiet",
            "--fastqout",
            null,
            "--fastqout_rev",
            null,
            "--fastqout_orphans",
            null,
            "--fastqout_orphans_rev",
            null,
        ];
        group.bench_function("vsearch", |b| b.iter(|| run(&vsearch, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
