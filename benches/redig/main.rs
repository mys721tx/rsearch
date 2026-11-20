//! Benchmarks of cutting sequences at a restriction pattern.
//!
//! `binary` runs the redig program end to end on a FASTA file of random
//! sequences of 1,000 bases. Most hold an EcoRI site (G^AATT_C) planted a
//! few times. It writes the fragments of both strands and the sequences
//! never cut. Set `VSEARCH` to a vsearch executable to time `vsearch --cut`
//! on the same file for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench redig
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// Sequences of `LENGTH` nucleotides.
const SEQUENCES: usize = 50_000;
const LENGTH: usize = 1_000;

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

/// Generates random sequences with up to three EcoRI sites each.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..SEQUENCES {
        let mut seq: Vec<u8> = (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect();
        for _ in 0..rng.below(4) {
            let at = rng.below(LENGTH - 6);
            seq[at..at + 6].copy_from_slice(b"GAATTC");
        }
        out.push_str(&format!(">seq{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    out
}

fn input_file() -> PathBuf {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("redig-bench.fa");
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
    let fa = input.to_str().unwrap();
    let redig = env!("CARGO_BIN_EXE_redig");
    let null = "/dev/null";

    let mut group = c.benchmark_group("redig/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let args = [
        "-i",
        fa,
        "--pattern",
        "G^AATT_C",
        "-o",
        null,
        "--rev-out",
        null,
        "--uncut",
        null,
    ];
    group.bench_function("redig", |b| b.iter(|| run(redig, &args)));
    if let Ok(vsearch) = std::env::var("VSEARCH") {
        let args = [
            "--cut",
            fa,
            "--cut_pattern",
            "G^AATT_C",
            "--quiet",
            "--fasta_width",
            "0",
            "--fastaout",
            null,
            "--fastaout_rev",
            null,
            "--fastaout_discarded",
            null,
        ];
        translated::bench(&mut group, "vsearch", &args);
        group.bench_function("vsearch", |b| b.iter(|| run(&vsearch, &args)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
