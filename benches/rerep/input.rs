//! The input of the rerep benchmarks, shared with compare/benches/rerep.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

use std::path::{Path, PathBuf};

/// Dereplicated reads of `LENGTH` nucleotides, with abundances from 1 to
/// `MAX_SIZE`.
pub const READS: usize = 20_000;
pub const LENGTH: usize = 250;
pub const MAX_SIZE: usize = 20;

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

/// Generates dereplicated reads as FASTQ and as FASTA, and returns them
/// with the number of copies they make.
pub fn generate() -> (String, String, usize) {
    let mut rng = Rng(0x5eed);
    let (mut fastq, mut fasta, mut copies) = (String::new(), String::new(), 0);
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        let size = 1 + rng.below(MAX_SIZE);
        copies += size;
        fastq.push_str(&format!("@read{i};size={size}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i};size={size}\n{seq}\n"));
    }
    (fastq, fasta, copies)
}

/// The generated input files, written once into the target directory, and
/// the number of copies they make.
pub struct Inputs {
    pub fastq: PathBuf,
    pub fasta: PathBuf,
    pub copies: usize,
}

pub fn input_files() -> Inputs {
    if let Some(fastq) = data::file("rerep", "uniques_fq") {
        let text = data::read(&fastq);
        let fasta = data::write("rerep", "uniques.fa", data::to_fasta(&data::parse(&text)));
        let copies = data::total_size(&text) as usize;
        return Inputs {
            fastq,
            fasta,
            copies,
        };
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fastq, fasta, copies) = generate();
    let inputs = Inputs {
        fastq: dir.join("rerep-bench.fq"),
        fasta: dir.join("rerep-bench.fa"),
        copies,
    };
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.fasta, fasta).unwrap();
    inputs
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("rerep", &["uniques_fq"])
}
