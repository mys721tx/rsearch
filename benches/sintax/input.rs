//! The input of the sintax benchmarks, shared with compare/benches/sintax.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub const READS: usize = 10_000;
pub const READ_LENGTH: usize = 250;
pub const DATABASE: usize = 5_000;
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

pub fn mutate(rng: &mut Rng, seq: &mut [u8], one_in: usize) {
    for base in seq.iter_mut() {
        if rng.below(one_in) == 0 {
            *base = b"ACGT"[rng.below(4)];
        }
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

/// Writes the reads, the FASTA database and its UDB file, built with
/// `mkudb`.
pub fn input_files(mkudb: &str) -> (PathBuf, PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut db = String::new();
    for i in 0..DATABASE {
        let t = rng.below(TEMPLATES);
        let mut seq = templates[t].clone();
        mutate(&mut rng, &mut seq, 20);
        db.push_str(&format!(
            ">d{i};tax=d:D{},p:P{},g:G{t}\n{}\n",
            t % 2,
            t % 7,
            String::from_utf8(seq).unwrap()
        ));
    }
    let mut reads = String::new();
    for i in 0..READS {
        let template = &templates[rng.below(TEMPLATES)];
        let start = rng.below(LENGTH - READ_LENGTH);
        let mut seq = template[start..start + READ_LENGTH].to_vec();
        mutate(&mut rng, &mut seq, 50);
        if rng.below(2) == 0 {
            seq = reverse_complement(&seq);
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    let (fa, dbfa, udb) = (
        tmp("sintax-bench.fa"),
        tmp("sintax-bench-db.fa"),
        tmp("sintax-bench-db.udb"),
    );
    std::fs::write(&fa, reads).unwrap();
    std::fs::write(&dbfa, db).unwrap();
    run(
        mkudb,
        &["-i", dbfa.to_str().unwrap(), "-o", udb.to_str().unwrap()],
    );
    (fa, dbfa, udb)
}

pub fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}
