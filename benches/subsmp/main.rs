//! Benchmarks of subsampling.
//!
//! `draw` measures `rsearch::random::draw` of half a million reads from a
//! million; `binary` runs the subsmp program end to end on a FASTQ file of
//! dereplicated reads: a tenth of the records, and a tenth of the reads by
//! abundance. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --fastx_subsample` on the same file for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench subsmp
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::random::{Mt19937_64, draw};

/// Records of `LENGTH` nucleotides, with abundances from 1 to 20.
const RECORDS: usize = 100_000;
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

/// Generates FASTQ records with abundances from 1 to 20.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut fastq = String::new();
    for i in 0..RECORDS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        let size = 1 + rng.below(20);
        fastq.push_str(&format!("@read{i};size={size}\n{seq}\n+\n{qual}\n"));
    }
    fastq
}

fn bench_draw(c: &mut Criterion) {
    let mut group = c.benchmark_group("subsmp/draw");
    let masses = vec![1u64; 1_000_000];
    group.throughput(Throughput::Elements(masses.len() as u64));
    group.bench_function("half_of_a_million", |b| {
        b.iter(|| black_box(draw(&masses, 500_000, &mut Mt19937_64::new(42))))
    });
    group.finish();
}

fn input_file() -> PathBuf {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("subsmp-bench.fq");
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
    let subsmp = env!("CARGO_BIN_EXE_subsmp");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("subsmp/binary");
    group.throughput(Throughput::Elements(RECORDS as u64));
    group.sample_size(10);
    for (name, by_size, vsearch_size) in [
        ("records", &[][..], &[][..]),
        ("by_size", &["--by-size"], &["--sizein", "--sizeout"]),
    ] {
        let args = [
            &[
                "-i",
                fq,
                "-o",
                "/dev/null",
                "--discarded",
                "/dev/null",
                "--seed",
                "42",
                "--pct",
                "10",
            ][..],
            by_size,
        ]
        .concat();
        group.bench_function(format!("subsmp/{name}"), |b| b.iter(|| run(subsmp, &args)));
        if let Some(vsearch) = &vsearch {
            let args = [
                &[
                    "--fastx_subsample",
                    fq,
                    "--quiet",
                    "--notrunclabels",
                    "--minseqlength",
                    "0",
                    "--fastqout",
                    "/dev/null",
                    "--fastqout_discarded",
                    "/dev/null",
                    "--randseed",
                    "42",
                    "--sample_pct",
                    "10",
                ][..],
                vsearch_size,
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

criterion_group!(benches, bench_draw, bench_binary);
criterion_main!(benches);
