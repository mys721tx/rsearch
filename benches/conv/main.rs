//! Benchmarks of quality conversion.
//!
//! `converter` measures `Converter::convert` on parsed reads; `binary` runs
//! the conv program end to end on files: Phred+33 to Phred+64 and back,
//! Solexa scores, limited output scores, relabelling and FASTA to FASTQ.
//! Set `VSEARCH` to a vsearch executable to time `vsearch --fastq_convert`
//! and `--fasta2fastq` on the same files for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench conv
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::quality::Conversion;

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

/// Generates Illumina-like reads with scores drifting from high to low
/// along the read, written with `offset` and scores from `lowest`.
fn generate(offset: i64, lowest: i64) -> String {
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
fn to_fasta(fastq: &str) -> String {
    fastq
        .lines()
        .collect::<Vec<_>>()
        .chunks(4)
        .map(|r| format!(">{}\n{}\n", &r[0][1..], r[1]))
        .collect()
}

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn bench_converter(c: &mut Criterion) {
    let reads = parse(&generate(33, 0));
    let mut group = c.benchmark_group("conv/converter");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, conversion) in [
        (
            "phred64",
            Conversion {
                offset_out: 64,
                max_out: 62,
                ..Conversion::default()
            },
        ),
        (
            "limited",
            Conversion {
                min_out: 2,
                max_out: 40,
                ..Conversion::default()
            },
        ),
    ] {
        let converter = conversion.converter();
        let mut out = Vec::new();
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    converter.convert(r.qual(), &mut out).unwrap();
                    black_box(&out);
                }
            })
        });
    }
    group.finish();
}

/// The generated input files, written once into the target directory.
struct Inputs {
    phred33: PathBuf,
    phred64: PathBuf,
    solexa: PathBuf,
    fasta: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        phred33: dir.join("conv-bench.fq"),
        phred64: dir.join("conv-bench-64.fq"),
        solexa: dir.join("conv-bench-solexa.fq"),
        fasta: dir.join("conv-bench.fa"),
    };
    let phred33 = generate(33, 0);
    std::fs::write(&inputs.fasta, to_fasta(&phred33)).unwrap();
    std::fs::write(&inputs.phred33, phred33).unwrap();
    std::fs::write(&inputs.phred64, generate(64, 0)).unwrap();
    std::fs::write(&inputs.solexa, generate(64, -5)).unwrap();
    inputs
}

/// A benchmark case: conv's arguments and the equivalent vsearch ones.
struct Case<'a> {
    name: &'static str,
    conv: Vec<&'a str>,
    vsearch: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.phred33.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let solexa = inputs.solexa.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let case = |name, file, conv: &[&'static str], command, vsearch: &[&'static str]| Case {
        name,
        conv: [&["-i", file, "-o", "/dev/null"], conv].concat(),
        vsearch: [
            &[command, file, "--fastqout", "/dev/null", "--quiet"],
            vsearch,
        ]
        .concat(),
    };
    let convert = "--fastq_convert";
    vec![
        case(
            "to_phred64",
            fq,
            &["--ascii-out", "64", "--qmax-out", "62"],
            convert,
            &["--fastq_asciiout", "64", "--fastq_qmaxout", "62"],
        ),
        case(
            "from_phred64",
            fq64,
            &["--ascii", "64"],
            convert,
            &["--fastq_ascii", "64"],
        ),
        case(
            "solexa",
            solexa,
            &["--ascii", "64", "--solexa"],
            convert,
            &["--fastq_ascii", "64", "--fastq_solexa"],
        ),
        case(
            "limited",
            fq,
            &["--qmin-out", "2", "--qmax-out", "40"],
            convert,
            &["--fastq_qminout", "2", "--fastq_qmaxout", "40"],
        ),
        case(
            "relabel",
            fq,
            &["--relabel", "R", "--sizeout"],
            convert,
            &["--relabel", "R", "--sizeout"],
        ),
        case("fasta", fa, &[], "--fasta2fastq", &[]),
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let conv = env!("CARGO_BIN_EXE_conv");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("conv/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("conv/{}", case.name), |b| {
            b.iter(|| run(conv, &case.conv))
        });
        if let Some(vsearch) = &vsearch {
            group.bench_function(format!("vsearch/{}", case.name), |b| {
                b.iter(|| run(vsearch, &case.vsearch))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_converter, bench_binary);
criterion_main!(benches);
