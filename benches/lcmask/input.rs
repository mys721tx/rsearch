//! The input of the lcmask benchmarks, shared with compare/benches/lcmask.

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

/// Generates FASTQ reads, a third with a dinucleotide repeat and a third
/// with a homopolymer of 20 somewhere in them.
pub fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..READS {
        let mut seq: Vec<u8> = (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect();
        let at = rng.below(LENGTH - 20);
        match i % 3 {
            0 => {
                let unit = [b"ACGT"[rng.below(4)], b"ACGT"[rng.below(4)]];
                for j in 0..20 {
                    seq[at + j] = unit[j % 2];
                }
            }
            1 => seq[at..at + 20].fill(b"ACGT"[rng.below(4)]),
            _ => {}
        }
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        let seq = String::from_utf8(seq).unwrap();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
    }
    fastq
}
pub fn input_file() -> PathBuf {
    if let Some(path) = data::file("lcmask", "reads") {
        return path;
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("lcmask-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("lcmask", &["reads"])
}

/// The reads as text: the dataset's, or the generated ones.
pub fn text() -> String {
    data::text("lcmask", "reads").unwrap_or_else(generate)
}
