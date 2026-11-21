//! The input of the udbstat benchmarks, shared with compare/benches/udbstat.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub const SEQUENCES: usize = 20_000;
pub const LENGTH: usize = 1_400;

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
}

pub fn tmp(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

/// Writes the FASTA input and builds its UDB file with `mkudb`.
pub fn input_file(mkudb: &str) -> PathBuf {
    let mut rng = Rng(0x5eed);
    let mut fasta = String::new();
    for i in 0..SEQUENCES {
        let seq: String = (0..LENGTH)
            .map(|_| b"ACGT"[(rng.next() % 4) as usize] as char)
            .collect();
        fasta.push_str(&format!(">seq{i}\n{seq}\n"));
    }
    let (fa, udb) = (tmp("udbstat-bench.fa"), tmp("udbstat-bench.udb"));
    std::fs::write(&fa, fasta).unwrap();
    let status = Command::new(mkudb)
        .args([
            "--word-len",
            "12",
            "-i",
            fa.to_str().unwrap(),
            "-o",
            udb.to_str().unwrap(),
        ])
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success());
    udb
}
