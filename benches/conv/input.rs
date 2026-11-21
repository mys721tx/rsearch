//! The input of the conv benchmarks, shared with compare/benches/conv.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::Command;

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

/// Generates Illumina-like reads with scores drifting from high to low
/// along the read, written with `offset` and scores from `lowest`.
pub fn generate(offset: i64, lowest: i64) -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        for _ in 0..LENGTH {
            q = (q + rng.below(4) as i64 - 2).clamp(lowest, 41);
            seq.push(char::from(b"ACGT"[rng.below(4)]));
            qual.push(char::from((q + offset) as u8));
        }
        out.push_str(&format!(
            "@read{i};size={}\n{seq}\n+\n{qual}\n",
            1 + rng.below(5)
        ));
    }
    out
}

/// Converts FASTQ text to FASTA.
pub fn to_fasta(fastq: &str) -> String {
    fastq
        .lines()
        .collect::<Vec<_>>()
        .chunks(4)
        .map(|r| format!(">{}\n{}\n", &r[0][1..], r[1]))
        .collect()
}

/// Pads `bytes` with zeros to a multiple of 8, as SFF sections are.
pub fn pad8(bytes: &mut Vec<u8>) {
    bytes.resize(bytes.len().div_ceil(8) * 8, 0);
}

/// Builds an SFF file of the reads of `fastq` (Phred+33), with 400 flows
/// per read, clip points 5 bases in from each end and an index at the end.
pub fn to_sff(fastq: &str) -> Vec<u8> {
    pub const FLOWS: u16 = 400;
    let header_length = (31 + FLOWS + 4).div_ceil(8) * 8;
    let lines: Vec<&str> = fastq.lines().collect();
    let mut body = Vec::new();
    for read in lines.chunks(4) {
        let (name, bases, quals) = (&read[0][1..], read[1].as_bytes(), read[3].as_bytes());
        let n = bases.len() as u16;
        let mut header = Vec::new();
        header.extend(
            (16 + name.len() as u16)
                .div_ceil(8)
                .wrapping_mul(8)
                .to_be_bytes(),
        );
        header.extend((name.len() as u16).to_be_bytes());
        header.extend(u32::from(n).to_be_bytes());
        for clip in [5, n - 5, 0, 0] {
            header.extend(clip.to_be_bytes());
        }
        header.extend(name.as_bytes());
        pad8(&mut header);
        body.extend(header);
        let mut data = vec![0; 2 * usize::from(FLOWS)];
        data.extend(vec![1; bases.len()]);
        data.extend(bases);
        data.extend(quals.iter().map(|q| q - 33));
        pad8(&mut data);
        body.extend(data);
    }
    let mut index = b".mft1.00".to_vec();
    index.extend([0; 8]);
    let index_offset = u64::from(header_length) + body.len() as u64;
    let mut header = b".sff".to_vec();
    header.extend(1u32.to_be_bytes());
    header.extend(index_offset.to_be_bytes());
    header.extend((index.len() as u32).to_be_bytes());
    header.extend((lines.len() as u32 / 4).to_be_bytes());
    header.extend(header_length.to_be_bytes());
    header.extend(4u16.to_be_bytes());
    header.extend(FLOWS.to_be_bytes());
    header.push(1);
    header.extend(b"TACG".repeat(usize::from(FLOWS) / 4));
    header.extend(b"TCAG");
    pad8(&mut header);
    pad8(&mut index);
    [header, body, index].concat()
}
/// The generated input files, written once into the target directory.
pub struct Inputs {
    pub phred33: PathBuf,
    pub phred64: PathBuf,
    pub solexa: PathBuf,
    pub fasta: PathBuf,
    pub sff: PathBuf,
    pub udb: PathBuf,
}

/// Writes the input files, the UDB file built with `mkudb`.
pub fn input_files(mkudb: &str) -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        phred33: dir.join("conv-bench.fq"),
        phred64: dir.join("conv-bench-64.fq"),
        solexa: dir.join("conv-bench-solexa.fq"),
        fasta: dir.join("conv-bench.fa"),
        sff: dir.join("conv-bench.sff"),
        udb: dir.join("conv-bench.udb"),
    };
    let phred33 = generate(33, 0);
    std::fs::write(&inputs.sff, to_sff(&phred33)).unwrap();
    std::fs::write(&inputs.fasta, to_fasta(&phred33)).unwrap();
    std::fs::write(&inputs.phred33, phred33).unwrap();
    std::fs::write(&inputs.phred64, generate(64, 0)).unwrap();
    std::fs::write(&inputs.solexa, generate(64, -5)).unwrap();
    let status = Command::new(mkudb)
        .args(["-i", inputs.fasta.to_str().unwrap()])
        .args(["-o", inputs.udb.to_str().unwrap()])
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success());
    inputs
}
