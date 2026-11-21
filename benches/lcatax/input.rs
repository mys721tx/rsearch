//! The input of the lcatax benchmarks, shared with compare/benches/lcatax.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

use std::path::{Path, PathBuf};

pub const READS: usize = 10_000;
pub const DISTINCT: usize = 20_000;
pub const DUPLICATES: usize = 5_000;
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

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

pub fn tmp(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

/// writes the reads and the database.
pub fn input_files() -> (PathBuf, PathBuf) {
    // Not the dataset's `db`. vsearch's command for the LCA aligns every
    // query with every target (`--id 1 --maxrejects 0`), hours on a
    // full-length reference. `lca_db`, a small reference with `;tax=`, is
    // searched against its own sequences. Identical sequences of different
    // taxa give the last common ancestor its work.
    if let Some(db) = data::file("lcatax", "lca_db") {
        let queries: Vec<data::Record> = data::parse(&data::read(&db))
            .into_iter()
            .enumerate()
            .map(|(i, r)| data::Record {
                label: format!("q{i}"),
                ..r
            })
            .collect();
        return (
            data::write("lcatax", "queries.fa", data::to_fasta(&queries)),
            db,
        );
    }
    let mut rng = Rng(0x5eed);
    let seqs: Vec<Vec<u8>> = (0..DISTINCT)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut db = String::new();
    for (i, seq) in seqs.iter().enumerate() {
        db.push_str(&format!(
            ">d{i};tax=d:D{},p:P{},g:G{}\n{}\n",
            i % 2,
            i % 7,
            i % 2000,
            String::from_utf8_lossy(seq)
        ));
    }
    for i in 0..DUPLICATES {
        let seq = &seqs[rng.below(DISTINCT)];
        db.push_str(&format!(">dup{i}\n{}\n", String::from_utf8_lossy(seq)));
    }
    let mut reads = String::new();
    for i in 0..READS {
        let mut seq = seqs[rng.below(DISTINCT)].clone();
        if rng.below(10) == 0 {
            let at = rng.below(LENGTH);
            seq[at] = b"ACGT"[rng.below(4)];
        }
        if rng.below(4) == 0 {
            seq = reverse_complement(&seq);
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8_lossy(&seq)));
    }
    let (fa, dbfa) = (tmp("lcatax-bench.fa"), tmp("lcatax-bench-db.fa"));
    std::fs::write(&fa, reads).unwrap();
    std::fs::write(&dbfa, db).unwrap();
    (fa, dbfa)
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("lcatax", &["lca_db"])
}
