//! The input of the mkudb benchmarks, shared with compare/benches/mkudb.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Sequences of `LENGTH` nucleotides, copies of `TEMPLATES` templates.
pub const SEQUENCES: usize = 20_000;
pub const LENGTH: usize = 1_400;
pub const TEMPLATES: usize = 50;

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

/// Generates copies of random templates with one base in twenty changed.
pub fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut out = String::new();
    for i in 0..SEQUENCES {
        let mut seq = templates[rng.below(TEMPLATES)].clone();
        for base in seq.iter_mut() {
            if rng.below(20) == 0 {
                *base = b"ACGT"[rng.below(4)];
            }
        }
        out.push_str(&format!(
            ">seq{i} ref\n{}\n",
            String::from_utf8(seq).unwrap()
        ));
    }
    out
}

pub fn tmp(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

/// Writes the FASTA input into the target directory.
pub fn input_file() -> PathBuf {
    let path = tmp("mkudb-bench.fa");
    std::fs::write(&path, generate()).unwrap();
    path
}
