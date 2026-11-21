//! The input of the dbsrch benchmarks, shared with compare/benches/dbsrch.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

pub const QUERIES: usize = 2_000;
pub const TARGETS: usize = 10_000;
pub const LENGTH: usize = 250;
pub const TEMPLATES: usize = 200;

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

/// the queries and the database.
pub fn input_files() -> (PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut fasta = |prefix: &str, n: usize| {
        let mut text = String::new();
        for i in 0..n {
            let mut seq = Vec::new();
            for &base in &templates[rng.below(TEMPLATES)] {
                match rng.below(60) {
                    0 => seq.push(b"ACGT"[rng.below(4)]),
                    1 => {}
                    2 => seq.extend([base, b"ACGT"[rng.below(4)]]),
                    _ => seq.push(base),
                }
            }
            text.push_str(&format!(
                ">{prefix}{i}\n{}\n",
                String::from_utf8(seq).unwrap()
            ));
        }
        text
    };
    let db = fasta("t", TARGETS);
    let q = fasta("q", QUERIES);
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (qp, dbp) = (
        dir.join("dbsrch-bench-q.fa"),
        dir.join("dbsrch-bench-db.fa"),
    );
    std::fs::write(&qp, q).unwrap();
    std::fs::write(&dbp, db).unwrap();
    (qp, dbp)
}
