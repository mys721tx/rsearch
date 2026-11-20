//! Benchmarks of DUST masking.
//!
//! `dust` measures `rsearch::dust::dust` on parsed reads, soft and hard;
//! `binary` runs the lcmask program end to end on a FASTQ file, soft and
//! hard. The reads are random with a repeat, a homopolymer or nothing
//! inserted, so that some windows are masked and most are searched in
//! full. Set `VSEARCH` to a vsearch executable to time `vsearch
//! --fastx_mask` on the same file for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench lcmask
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::dust::{Style, dust};
use rsearch::fastx::{Reader, Record};

/// Reads of `LENGTH` nucleotides.
const READS: usize = 100_000;
const LENGTH: usize = 250;

/// A small deterministic generator (xorshift64*). Every run benchmarks the
/// same data without extra dependencies.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// Generates FASTQ reads, a third with a dinucleotide repeat and a third
/// with a homopolymer of 20 somewhere in them.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..READS {
        let mut seq: Vec<u8> = (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect();
        let at = rng.below(LENGTH - 20);
        match i % 3 {
            0 => {
                let unit = [b"ACGT"[rng.below(4)], b"ACGT"[rng.below(4)]];
                for j in 0..20 {
                    seq[at + j] = unit[j % 2];
                }
            }
            1 => seq[at..at + 20].fill(b"ACGT"[rng.below(4)]),
            _ => {}
        }
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        let seq = String::from_utf8(seq).unwrap();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
    }
    fastq
}

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_dust(c: &mut Criterion) {
    let reads = parse(&generate());
    let mut group = c.benchmark_group("lcmask/dust");
    group.throughput(Throughput::Elements(READS as u64));
    let mut seq = Vec::new();
    for (name, style) in [("soft", Style::Soft), ("hard", Style::Hard)] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    seq.clear();
                    seq.extend_from_slice(r.seq());
                    dust(&mut seq, style);
                    black_box(&seq);
                }
            })
        });
    }
    group.finish();
}

fn input_file() -> PathBuf {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("lcmask-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fq = input.to_str().unwrap();
    let lcmask = env!("CARGO_BIN_EXE_lcmask");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("lcmask/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, hard, hardmask) in [
        ("soft", &[][..], &[][..]),
        ("hard", &["--hard"], &["--hardmask"]),
    ] {
        let args = [&["-i", fq, "-o", "/dev/null"][..], hard].concat();
        group.bench_function(format!("lcmask/{name}"), |b| b.iter(|| run(lcmask, &args)));
        if let Some(vsearch) = &vsearch {
            let args = [
                &[
                    "--fastx_mask",
                    fq,
                    "--quiet",
                    "--notrunclabels",
                    "--minseqlength",
                    "0",
                    "--threads",
                    "1",
                    "--fastqout",
                    "/dev/null",
                ][..],
                hardmask,
            ]
            .concat();
            translated::bench(&mut group, format!("vsearch/{name}"), &args);
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_dust, bench_binary);
criterion_main!(benches);
