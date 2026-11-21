//! The input of the allpr benchmarks, shared with compare/benches/allpr.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

pub const SEQUENCES: usize = 600;
pub const LENGTH: usize = 250;
pub const TEMPLATES: usize = 20;

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

pub fn input_file() -> PathBuf {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut fasta = String::new();
    for i in 0..SEQUENCES {
        let mut seq = Vec::new();
        for &base in &templates[rng.below(TEMPLATES)] {
            match rng.below(60) {
                0 => seq.push(b"ACGT"[rng.below(4)]),
                1 => {}
                2 => seq.extend([base, b"ACGT"[rng.below(4)]]),
                _ => seq.push(base),
            }
        }
        fasta.push_str(&format!(">s{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("allpr-bench.fa");
    std::fs::write(&path, fasta).unwrap();
    path
}
