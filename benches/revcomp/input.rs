//! The input of the revcomp benchmarks, shared with compare/benches/revcomp.

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

/// Generates reads as FASTQ and as FASTA, mostly ACGT with some lower case
/// and ambiguity codes.
pub fn generate() -> (String, String) {
    pub const CODES: &[u8] = b"ACGTACGTACGTACGTACGTACGTacgtNRY";
    let mut rng = Rng(0x5eed);
    let (mut fastq, mut fasta) = (String::new(), String::new());
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(CODES[rng.below(CODES.len())]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i}\n{seq}\n"));
    }
    (fastq, fasta)
}
/// The generated input files, written once into the target directory.
pub struct Inputs {
    pub fastq: PathBuf,
    pub fasta: PathBuf,
}

pub fn input_files() -> Inputs {
    if let Some(fastq) = data::file("revcomp", "reads") {
        let reads = data::parse(&data::read(&fastq));
        let fasta = data::write("revcomp", "reads.fa", data::to_fasta(&reads));
        return Inputs { fastq, fasta };
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fastq, fasta) = generate();
    let inputs = Inputs {
        fastq: dir.join("revcomp-bench.fq"),
        fasta: dir.join("revcomp-bench.fa"),
    };
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.fasta, fasta).unwrap();
    inputs
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("revcomp", &["reads"])
}

/// The reads as FASTQ text: the dataset's, or the generated ones.
pub fn text() -> String {
    data::text("revcomp", "reads").unwrap_or_else(|| generate().0)
}
