//! The input of the uchime benchmarks, shared with compare/benches/uchime.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

pub const SEQUENCES: usize = 4_000;
pub const LENGTH: usize = 250;
pub const PARENTS: usize = 100;

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

/// the input, sorted by decreasing abundance, the queries (the input but
/// the parents) and the parents.
pub fn input_files(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let parents: Vec<Vec<u8>> = (0..PARENTS)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let (mut all, mut queries, mut db) = (String::new(), String::new(), String::new());
    for (i, p) in parents.iter().enumerate() {
        let record = format!(
            ">p{i};size={}\n{}\n",
            100_000 - i,
            String::from_utf8_lossy(p)
        );
        all.push_str(&record);
        db.push_str(&record);
    }
    for i in PARENTS..SEQUENCES {
        let a = rng.below(PARENTS);
        let mut seq = parents[a].clone();
        // a third of them chimeras
        if rng.below(3) == 0 {
            let b = rng.below(PARENTS);
            let cut = 50 + rng.below(LENGTH - 100);
            seq[cut..].copy_from_slice(&parents[b][cut..]);
        }
        for base in seq.iter_mut() {
            if rng.below(100) == 0 {
                *base = b"ACGT"[rng.below(4)];
            }
        }
        // decreasing abundances
        let record = format!(
            ">s{i};size={}\n{}\n",
            50_000 - i,
            String::from_utf8_lossy(&seq)
        );
        all.push_str(&record);
        queries.push_str(&record);
    }
    let paths = ["all", "queries", "db"].map(|n| dir.join(format!("uchime-bench-{n}.fa")));
    for (path, text) in paths.iter().zip([&all, &queries, &db]) {
        std::fs::write(path, text).unwrap();
    }
    let [a, q, d] = paths;
    (a, q, d)
}
