//! Benchmarks of read trimming.
//!
//! `trim` measures `Filter::apply` with trimming settings on parsed reads;
//! `binary` runs the sqtrim program end to end on files: stripping ends,
//! truncating to a length, at a low score, at an expected-error limit or
//! rate, all of them together, Phred+64 and FASTA. Set `VSEARCH` to a
//! vsearch executable to time `vsearch --fastx_filter` with the same
//! trimming options on the same files for comparison, writing its discarded
//! reads too, since sqtrim writes every read:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench sqtrim
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::filter::Filter;

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
/// along the read, written with `offset`.
fn generate(offset: i64) -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        for _ in 0..LENGTH {
            q = (q + rng.below(4) as i64 - 2).clamp(2, 41);
            seq.push(char::from(b"ACGT"[rng.below(4)]));
            qual.push(char::from((q + offset) as u8));
        }
        out.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
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

fn bench_trim(c: &mut Criterion) {
    let reads = parse(&generate(33));
    let mut group = c.benchmark_group("sqtrim/trim");
    group.throughput(Throughput::Elements(READS as u64));
    let trimming = Filter {
        min_len: 0,
        ..Filter::default()
    };
    for (name, filter) in [
        (
            "strip_len",
            Filter {
                strip_left: 10,
                trunc_keep: Some(200),
                ..trimming.clone()
            },
        ),
        (
            "qual",
            Filter {
                trunc_qual: Some(10),
                ..trimming.clone()
            },
        ),
        (
            "max_ee",
            Filter {
                trunc_ee: Some(1.0),
                ..trimming.clone()
            },
        ),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for r in &reads {
                    black_box(filter.apply(r.seq(), Some(r.qual()), 1).unwrap());
                }
            })
        });
    }
    group.finish();
}

/// The generated input files, written once into the target directory.
struct Inputs {
    fastq: PathBuf,
    phred64: PathBuf,
    fasta: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fastq: dir.join("sqtrim-bench.fq"),
        phred64: dir.join("sqtrim-bench-64.fq"),
        fasta: dir.join("sqtrim-bench.fa"),
    };
    let fastq = generate(33);
    std::fs::write(&inputs.fasta, to_fasta(&fastq)).unwrap();
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.phred64, generate(64)).unwrap();
    inputs
}

/// A benchmark case: sqtrim's arguments and the equivalent vsearch ones.
struct Case<'a> {
    name: &'static str,
    sqtrim: Vec<&'a str>,
    vsearch: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let case = |name, file, sqtrim: &[&'static str], vsearch: &[&'static str]| {
        let output: &[&str] = if file == fa {
            &[
                "--fastaout",
                "/dev/null",
                "--fastaout_discarded",
                "/dev/null",
                "--fasta_width",
                "0",
            ]
        } else {
            &[
                "--fastqout",
                "/dev/null",
                "--fastqout_discarded",
                "/dev/null",
            ]
        };
        Case {
            name,
            sqtrim: [&["-i", file, "-o", "/dev/null"], sqtrim].concat(),
            vsearch: [&["--fastx_filter", file, "--quiet"], output, vsearch].concat(),
        }
    };
    vec![
        case(
            "strip",
            fq,
            &["--strip-left", "10", "--strip-right", "10"],
            &["--fastq_stripleft", "10", "--fastq_stripright", "10"],
        ),
        case(
            "len",
            fq,
            &["--len", "200"],
            &["--fastq_trunclen_keep", "200"],
        ),
        case("qual", fq, &["--qual", "10"], &["--fastq_truncqual", "10"]),
        case("max_ee", fq, &["--max-ee", "1"], &["--fastq_truncee", "1"]),
        case(
            "max_ee_rate",
            fq,
            &["--max-ee-rate", "0.01"],
            &["--fastq_truncee_rate", "0.01"],
        ),
        case(
            "all",
            fq,
            &[
                "--strip-left",
                "5",
                "--len",
                "200",
                "--qual",
                "5",
                "--max-ee",
                "2",
            ],
            &[
                "--fastq_stripleft",
                "5",
                "--fastq_trunclen_keep",
                "200",
                "--fastq_truncqual",
                "5",
                "--fastq_truncee",
                "2",
            ],
        ),
        case(
            "phred64",
            fq64,
            &["--ascii", "64", "--qual", "10"],
            &["--fastq_ascii", "64", "--fastq_truncqual", "10"],
        ),
        case(
            "fasta",
            fa,
            &["--strip-left", "10", "--len", "200"],
            &["--fastq_stripleft", "10", "--fastq_trunclen_keep", "200"],
        ),
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let sqtrim = env!("CARGO_BIN_EXE_sqtrim");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("sqtrim/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("sqtrim/{}", case.name), |b| {
            b.iter(|| run(sqtrim, &case.sqtrim))
        });
        if let Some(vsearch) = &vsearch {
            translated::bench(&mut group, format!("vsearch/{}", case.name), &case.vsearch);
            group.bench_function(format!("vsearch/{}", case.name), |b| {
                b.iter(|| run(vsearch, &case.vsearch))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_trim, bench_binary);
criterion_main!(benches);
