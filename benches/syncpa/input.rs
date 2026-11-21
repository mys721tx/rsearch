//! The input of the syncpa benchmarks, shared with compare/benches/syncpa.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// Pairs of reads of `LENGTH` nucleotides.
pub const PAIRS: usize = 100_000;
pub const LENGTH: usize = 150;

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

/// A FASTQ read with Illumina-like headers ending in /1 or /2.
pub fn read(rng: &mut Rng, i: usize, mate: u8) -> String {
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
pub fn generate() -> (String, String) {
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

pub fn input_files() -> (PathBuf, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (forward, reverse) = generate();
    let paths = (dir.join("syncpa-bench-1.fq"), dir.join("syncpa-bench-2.fq"));
    std::fs::write(&paths.0, forward).unwrap();
    std::fs::write(&paths.1, reverse).unwrap();
    paths
}
