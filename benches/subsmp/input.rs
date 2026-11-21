//! The input of the subsmp benchmarks, shared with compare/benches/subsmp.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Records of `LENGTH` nucleotides, with abundances from 1 to 20.
pub const RECORDS: usize = 100_000;
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

/// Generates FASTQ records with abundances from 1 to 20.
pub fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..RECORDS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        let size = 1 + rng.below(20);
        fastq.push_str(&format!("@read{i};size={size}\n{seq}\n+\n{qual}\n"));
    }
    fastq
}
pub fn input_file() -> PathBuf {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("subsmp-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}
