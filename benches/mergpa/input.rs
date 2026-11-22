//! The input of the mergpa benchmarks, shared with compare/benches/mergpa.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

use std::path::{Path, PathBuf};

/// Pairs of reads of `LENGTH` nucleotides.
pub const PAIRS: usize = 100_000;
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

/// A read's sequence and scores.
pub type Read = (Vec<u8>, Vec<u8>);

/// Scores around Q30, from Q2 to Q41, with an error now and then.
pub fn read(rng: &mut Rng, bases: impl Iterator<Item = u8>) -> Read {
    let mut seq = Vec::with_capacity(LENGTH);
    let mut qual = Vec::with_capacity(LENGTH);
    for base in bases {
        let score = (18 + rng.below(23)).max(2) as u8;
        let base = if rng.below(200) == 0 {
            b"ACGT"[rng.below(4)]
        } else {
            base
        };
        seq.push(base);
        qual.push(33 + score);
    }
    (seq, qual)
}

/// Generates the forward and reverse reads of each fragment.
pub fn generate() -> Vec<(Read, Read)> {
    let mut rng = Rng(0x5eed);
    (0..PAIRS)
        .map(|_| {
            let length = 300 + rng.below(151);
            let fragment: Vec<u8> = (0..length).map(|_| b"ACGT"[rng.below(4)]).collect();
            let forward = read(&mut rng, fragment[..LENGTH].iter().copied());
            let reverse = read(
                &mut rng,
                fragment[length - LENGTH..].iter().rev().map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                }),
            );
            (forward, reverse)
        })
        .collect()
}

/// Pairs of fragments that hold a tandem repeat of 40 bases or more, a unit
/// of 2 to 8 repeated, where both reads overlap. Many diagonals then share
/// enough 5-mers to be scored, and most fall apart along the way.
pub fn generate_repeats() -> Vec<(Read, Read)> {
    let mut rng = Rng(0x7e9ea7);
    (0..PAIRS / 4)
        .map(|_| {
            // an overlap of 60 to 200 bases, from length - 250 to 250
            let length = 300 + rng.below(141);
            let overlap = 2 * LENGTH - length;
            let unit: Vec<u8> = (0..2 + rng.below(7))
                .map(|_| b"ACGT"[rng.below(4)])
                .collect();
            // a repeat of 40 bases to 10 short of the overlap, within it
            let span = 40 + rng.below(overlap - 49);
            let start = length - LENGTH + rng.below(overlap - span + 1);
            let fragment: Vec<u8> = (0..length)
                .map(|i| {
                    if (start..start + span).contains(&i) {
                        unit[(i - start) % unit.len()]
                    } else {
                        b"ACGT"[rng.below(4)]
                    }
                })
                .collect();
            let forward = read(&mut rng, fragment[..LENGTH].iter().copied());
            let reverse = read(
                &mut rng,
                fragment[length - LENGTH..].iter().rev().map(|&b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                }),
            );
            (forward, reverse)
        })
        .collect()
}

/// The files of the repeat pairs, written out.
pub fn repeat_files() -> (PathBuf, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let paths = (
        dir.join("mergpa-bench-repeats-1.fq"),
        dir.join("mergpa-bench-repeats-2.fq"),
    );
    let pairs = generate_repeats();
    std::fs::write(
        &paths.0,
        fastq(pairs.iter().map(|p| p.0.clone()).enumerate(), 1),
    )
    .unwrap();
    std::fs::write(
        &paths.1,
        fastq(pairs.iter().map(|p| p.1.clone()).enumerate(), 2),
    )
    .unwrap();
    paths
}

pub fn fastq(reads: impl Iterator<Item = (usize, Read)>, mate: u8) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, (seq, qual)) in reads {
        out.extend_from_slice(format!("@read{i}/{mate}\n").as_bytes());
        out.extend_from_slice(&seq);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&qual);
        out.push(b'\n');
    }
    out
}

pub fn input_files(pairs: &[(Read, Read)]) -> (PathBuf, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let paths = (dir.join("mergpa-bench-1.fq"), dir.join("mergpa-bench-2.fq"));
    let forward = fastq(pairs.iter().map(|p| p.0.clone()).enumerate(), 1);
    let reverse = fastq(pairs.iter().map(|p| p.1.clone()).enumerate(), 2);
    std::fs::write(&paths.0, forward).unwrap();
    std::fs::write(&paths.1, reverse).unwrap();
    paths
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("mergpa", &["reads", "mate"])
}

/// The pairs to merge in memory: the dataset's, or the generated ones.
pub fn pairs() -> Vec<(Read, Read)> {
    match (data::file("mergpa", "reads"), data::file("mergpa", "mate")) {
        (Some(f), Some(r)) => {
            let read = |rec: data::Record| (rec.seq.into_bytes(), rec.qual.unwrap().into_bytes());
            let f = data::parse(&data::read(&f));
            let r = data::parse(&data::read(&r));
            f.into_iter()
                .map(read)
                .zip(r.into_iter().map(read))
                .collect()
        }
        _ => generate(),
    }
}

/// The files of pairs: the dataset's, or the generated ones written out.
pub fn pair_files() -> (PathBuf, PathBuf) {
    match (data::file("mergpa", "reads"), data::file("mergpa", "mate")) {
        (Some(f), Some(r)) => (f, r),
        _ => input_files(&generate()),
    }
}
