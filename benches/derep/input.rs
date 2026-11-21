//! The input of the derep benchmarks, shared with compare/benches/derep.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

use std::path::{Path, PathBuf};

/// Reads drawn from `UNIQUES` distinct sequences of `LENGTH` nucleotides.
pub const READS: usize = 100_000;
pub const UNIQUES: usize = 2_000;
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

/// The reverse complement of a sequence of upper case ACGT.
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

/// Generates reads with a skewed abundance distribution: read i comes from
/// a unique sequence chosen with a bias towards low indices. A tenth of the
/// reads are lower case, and a tenth are reverse complemented.
pub fn generate(fastq: bool) -> String {
    let mut rng = Rng(0x5eed);
    let uniques: Vec<Vec<u8>> = (0..UNIQUES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut out = String::new();
    for i in 0..READS {
        // the minimum of two uniform draws favours abundant sequences
        let u = rng.below(UNIQUES).min(rng.below(UNIQUES));
        let mut seq = uniques[u].clone();
        match rng.below(10) {
            0 => seq.make_ascii_lowercase(),
            1 => seq = reverse_complement(&seq),
            _ => {}
        }
        let seq = String::from_utf8(seq).unwrap();
        if fastq {
            let qual: String = (0..LENGTH)
                .map(|_| char::from(b'#' + rng.below(39) as u8))
                .collect();
            out.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        } else {
            out.push_str(&format!(">read{i}\n{seq}\n"));
        }
    }
    out
}
/// Re-encodes Phred+33 FASTQ text as Phred+64.
pub fn to_phred64(fastq: &str) -> String {
    fastq
        .lines()
        .enumerate()
        .map(|(i, line)| {
            if i % 4 == 3 {
                line.bytes().map(|q| char::from(q + 31)).collect::<String>() + "\n"
            } else {
                format!("{line}\n")
            }
        })
        .collect()
}

/// The generated input files: FASTA, FASTQ, and the FASTQ as Phred+64.
pub struct Inputs {
    pub fasta: PathBuf,
    pub fastq: PathBuf,
    pub phred64: PathBuf,
}

/// Writes the generated input files once, into the target directory.
pub fn input_files() -> Inputs {
    if let Some(fastq) = data::file("derep", "reads") {
        let text = data::read(&fastq);
        return Inputs {
            fasta: data::write("derep", "reads.fa", data::to_fasta(&data::parse(&text))),
            phred64: data::write("derep", "reads-64.fq", data::phred64(&text)),
            fastq,
        };
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fasta: dir.join("derep-bench.fa"),
        fastq: dir.join("derep-bench.fq"),
        phred64: dir.join("derep-bench-64.fq"),
    };
    let fastq = generate(true);
    std::fs::write(&inputs.fasta, generate(false)).unwrap();
    std::fs::write(&inputs.phred64, to_phred64(&fastq)).unwrap();
    std::fs::write(&inputs.fastq, fastq).unwrap();
    inputs
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("derep", &["reads"])
}

/// The reads as text, FASTQ or FASTA: the dataset's, or the generated
/// ones.
pub fn text(fastq: bool) -> String {
    match data::text("derep", "reads") {
        Some(text) if fastq => text,
        Some(text) => data::to_fasta(&data::parse(&text)),
        None => generate(fastq),
    }
}
