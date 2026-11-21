//! The input of the redig benchmarks, shared with compare/benches/redig.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

use std::path::{Path, PathBuf};

/// Sequences of `LENGTH` nucleotides.
pub const SEQUENCES: usize = 50_000;
pub const LENGTH: usize = 1_000;

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

/// Generates random sequences with up to three EcoRI sites each.
pub fn generate() -> String {
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

pub fn input_file() -> PathBuf {
    if let Some(path) = data::file("redig", "db") {
        return path;
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("redig-bench.fa");
    std::fs::write(&path, generate()).unwrap();
    path
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("redig", &["db"])
}
