//! The input of the scrmbl benchmarks, shared with compare/benches/scrmbl.

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

/// Generates FASTQ reads.
pub fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
    }
    fastq
}
pub fn input_file() -> PathBuf {
    if let Some(path) = data::file("scrmbl", "reads") {
        return path;
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("scrmbl-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("scrmbl", &["reads"])
}

/// The reads as text: the dataset's, or the generated ones.
pub fn text() -> String {
    data::text("scrmbl", "reads").unwrap_or_else(generate)
}
