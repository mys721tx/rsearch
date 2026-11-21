//! The input of the eegrid benchmarks, shared with compare/benches/eegrid.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

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

/// Generates Illumina-like reads with scores drifting from high to low
/// along the read. A read ends in a run of score 2 now and then. One read
/// in three is trimmed to a length from 150 to 250.
pub fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        let tail = if rng.below(10) == 0 { 20 } else { 0 };
        let length = if i % 3 == 0 {
            LENGTH - rng.below(101)
        } else {
            LENGTH
        };
        for j in 0..length {
            q = (q + rng.below(4) as i64 - 2).clamp(2, 41);
            seq.push(char::from(b"ACGTACGTACGTN"[rng.below(13)]));
            let score = if j + tail >= LENGTH { 2 } else { q };
            qual.push(char::from((score + 33) as u8));
        }
        out.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
    }
    out
}

pub fn input_file() -> PathBuf {
    if let Some(path) = data::file("eegrid", "reads") {
        return path;
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("eegrid-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("eegrid", &["reads"])
}
