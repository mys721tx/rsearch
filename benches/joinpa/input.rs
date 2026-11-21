//! The input of the joinpa benchmarks: forward FASTQ reads and their mates,
//! shared with compare/benches/joinpa.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Pairs of reads of `LENGTH` nucleotides.
pub const PAIRS: usize = 100_000;
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

/// A FASTQ read with an Illumina-like header ending in /1 or /2.
pub fn read(rng: &mut Rng, i: usize, mate: u8) -> String {
    let seq: String = (0..LENGTH)
        .map(|_| char::from(b"ACGTACGTACGTN"[rng.below(13)]))
        .collect();
    let qual: String = (0..LENGTH)
        .map(|_| char::from(35 + rng.below(39) as u8))
        .collect();
    format!("@M00123:45:000000000-ABCDE:1:1101:{i}:1000/{mate}\n{seq}\n+\n{qual}\n")
}

pub fn input_files() -> (PathBuf, PathBuf) {
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
