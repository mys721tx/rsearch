//! The input of the sortby benchmarks, shared with compare/benches/sortby.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Reads of 150 to 300 nucleotides.
pub const READS: usize = 100_000;

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

/// Generates reads as FASTQ and as FASTA, with lengths from 150 to 300 and
/// abundances from 1 to 1000, so that both keys have ties.
pub fn generate() -> (String, String) {
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
pub struct Inputs {
    pub fastq: PathBuf,
    pub fasta: PathBuf,
}

pub fn input_files() -> Inputs {
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
