//! The input of the unoise benchmarks, shared with compare/benches/unoise.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

use std::path::{Path, PathBuf};

pub const SEQUENCES: usize = 10_000;
pub const LENGTH: usize = 250;
pub const TEMPLATES: usize = 300;

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
    if let Some(path) = data::file("unoise", "uniques") {
        return path;
    }
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut fasta = String::new();
    for i in 0..SEQUENCES {
        let template = &templates[if i < TEMPLATES {
            i
        } else {
            rng.below(TEMPLATES)
        }];
        let mut seq = Vec::new();
        for &base in template {
            match if i < TEMPLATES { 3 } else { rng.below(300) } {
                0 => seq.push(b"ACGT"[rng.below(4)]),
                1 => {}
                2 => seq.extend([base, b"ACGT"[rng.below(4)]]),
                _ => seq.push(base),
            }
        }
        // decreasing abundances, the templates first
        let size = 2 * SEQUENCES - i;
        fasta.push_str(&format!(
            ">s{i};size={size}\n{}\n",
            String::from_utf8(seq).unwrap()
        ));
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("unoise-bench.fa");
    std::fs::write(&path, fasta).unwrap();
    path
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("unoise", &["uniques"])
}
