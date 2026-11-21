//! The input of the relab benchmarks, shared with compare/benches/relab.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Reads of `LENGTH` nucleotides.
pub const READS: usize = 100_000;
pub const LENGTH: usize = 250;

/// A small deterministic generator (xorshift64*). Every run benchmarks the
/// same data without extra dependencies.
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    pub fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Generates Illumina-like Phred+33 reads with scores drifting from high to
/// low along the read, and headers carrying the three annotations relab
/// removes.
pub fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        for _ in 0..LENGTH {
            q = (q + rng.below(4) as i64 - 2).clamp(2, 41);
            seq.push(char::from(b"ACGT"[rng.below(4)]));
            qual.push(char::from((q + 33) as u8));
        }
        out.push_str(&format!(
            "@read{i};size={};ee=0.5;length={LENGTH}\n{seq}\n+\n{qual}\n",
            1 + rng.below(5)
        ));
    }
    out
}

/// Converts FASTQ text to FASTA.
pub fn to_fasta(fastq: &str) -> String {
    fastq
        .lines()
        .collect::<Vec<_>>()
        .chunks(4)
        .map(|r| format!(">{}\n{}\n", &r[0][1..], r[1]))
        .collect()
}
/// The generated input files, written once into the target directory.
pub struct Inputs {
    pub fastq: PathBuf,
    pub fasta: PathBuf,
}

pub fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fastq: dir.join("relab-bench.fq"),
        fasta: dir.join("relab-bench.fa"),
    };
    let fastq = generate();
    std::fs::write(&inputs.fasta, to_fasta(&fastq)).unwrap();
    std::fs::write(&inputs.fastq, fastq).unwrap();
    inputs
}
