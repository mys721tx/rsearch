//! The input of the getseq benchmarks, shared with compare/benches/getseq.

#![allow(dead_code)]

#[path = "../common/data.rs"]
pub mod data;

use std::path::{Path, PathBuf};

/// Reads of `LENGTH` nucleotides, from `SAMPLES` samples.
pub const READS: usize = 100_000;
pub const LENGTH: usize = 250;
pub const SAMPLES: usize = 100;

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

/// The header of read `i`.
pub fn header(i: usize) -> String {
    format!("read{i};sample=s{};size=1", i % SAMPLES)
}

/// Generates FASTQ reads with sample annotations.
pub fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@{}\n{seq}\n+\n{qual}\n", header(i)));
    }
    fastq
}
/// The generated input files, written once into the target directory.
pub struct Inputs {
    pub fastq: PathBuf,
    pub labels: PathBuf,
    pub words: PathBuf,
}

pub fn input_files() -> Inputs {
    if let Some(fastq) = data::file("getseq", "reads") {
        let labels: String = real_labels().into_iter().map(|l| l + "\n").collect();
        // the first words of the labels (run accessions in ENA's labels),
        // every other one
        let mut words: Vec<String> = Vec::new();
        for r in data::parse(&data::read(&fastq)) {
            let word: String = r
                .label
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                .collect();
            if !words.contains(&word) {
                words.push(word);
            }
        }
        let words: String = words.iter().step_by(2).map(|w| format!("{w}\n")).collect();
        return Inputs {
            fastq,
            labels: data::write("getseq", "labels.txt", labels),
            words: data::write("getseq", "words.txt", words),
        };
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fastq: dir.join("getseq-bench.fq"),
        labels: dir.join("getseq-bench-labels.txt"),
        words: dir.join("getseq-bench-words.txt"),
    };
    std::fs::write(&inputs.fastq, generate()).unwrap();
    let labels: String = (0..READS).step_by(10).map(|i| header(i) + "\n").collect();
    std::fs::write(&inputs.labels, labels).unwrap();
    let words: String = (0..SAMPLES).step_by(7).map(|s| format!("s{s}\n")).collect();
    std::fs::write(&inputs.words, words).unwrap();
    inputs
}
/// writes reads drawn from 20,000 sequences, a tenth of them changed, and
/// a database of those sequences and 5,000 duplicates.
pub fn hit_files() -> (PathBuf, PathBuf) {
    if let Some(fastq) = data::file("getseq", "reads") {
        let reads = data::parse(&data::read(&fastq));
        return (
            data::write("getseq", "hits.fa", data::to_fasta(&reads)),
            data::write("getseq", "hits-db.fa", data::exact_db(&reads)),
        );
    }
    let mut rng = Rng(0x5eed);
    let seqs: Vec<String> = (0..20_000)
        .map(|_| (0..150).map(|_| b"ACGT"[rng.below(4)] as char).collect())
        .collect();
    let mut db = String::new();
    for (i, seq) in seqs.iter().enumerate() {
        db.push_str(&format!(">d{i}\n{seq}\n"));
    }
    for i in 0..5_000 {
        db.push_str(&format!(">dup{i}\n{}\n", seqs[rng.below(seqs.len())]));
    }
    let mut reads = String::new();
    for i in 0..100_000 {
        let mut seq = seqs[rng.below(seqs.len())].clone().into_bytes();
        if rng.below(10) == 0 {
            let at = rng.below(150);
            seq[at] = b"ACGT"[rng.below(4)];
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fa, dbfa) = (dir.join("getseq-hits.fa"), dir.join("getseq-hits-db.fa"));
    std::fs::write(&fa, reads).unwrap();
    std::fs::write(&dbfa, db).unwrap();
    (fa, dbfa)
}

/// Whether a dataset is set that lacks the files this bench reads.
pub fn skip() -> bool {
    data::skip("getseq", &["reads"])
}

/// Every tenth label of the dataset's reads.
pub fn real_labels() -> Vec<String> {
    let text = data::text("getseq", "reads").unwrap_or_default();
    data::parse(&text)
        .into_iter()
        .step_by(10)
        .map(|r| r.label)
        .collect()
}

/// The reads as text: the dataset's, or the generated ones.
pub fn text() -> String {
    data::text("getseq", "reads").unwrap_or_else(generate)
}

/// The labels the in-memory selection looks for: every tenth read's.
pub fn labels() -> Vec<String> {
    if data::id().is_some() {
        real_labels()
    } else {
        (0..READS).step_by(10).map(header).collect()
    }
}
