//! Benchmarks of label rewriting.
//!
//! `format` measures `Label::format` on parsed reads; `binary` runs the
//! relab program end to end on files: copying, relabelling by prefix,
//! sequence, SHA1 and MD5, removing and adding annotations, and expected
//! errors, on FASTQ and FASTA. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --fastx_filter` with the same label options on the same files
//! for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench relab
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::label::{Label, Relabel, abundance};

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

/// Generates Illumina-like Phred+33 reads with scores drifting from high to
/// low along the read, and headers carrying the three annotations relab
/// removes.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        for _ in 0..LENGTH {
            q = (q + rng.below(4) as i64 - 2).clamp(2, 41);
            seq.push(char::from(b"ACGT"[rng.below(4)]));
            qual.push(char::from((q + 33) as u8));
        }
        out.push_str(&format!(
            "@read{i};size={};ee=0.5;length={LENGTH}\n{seq}\n+\n{qual}\n",
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

fn bench_format(c: &mut Criterion) {
    let reads = parse(&generate());
    let mut group = c.benchmark_group("relab/format");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, label) in [
        (
            "prefix",
            Label {
                relabel: Relabel::Prefix("R".into()),
                sizeout: true,
                ..Label::default()
            },
        ),
        (
            "sha1",
            Label {
                relabel: Relabel::Sha1,
                ..Label::default()
            },
        ),
        (
            "strip",
            Label {
                xsize: true,
                xee: true,
                xlength: true,
                ..Label::default()
            },
        ),
        (
            "annotate",
            Label {
                sizeout: true,
                lengthout: true,
                sample: Some("S1".into()),
                ..Label::default()
            },
        ),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| {
                for (i, r) in reads.iter().enumerate() {
                    let size = abundance(r.header());
                    black_box(label.format(r.header(), r.seq(), size, i as u64 + 1, None));
                }
            })
        });
    }
    group.finish();
}

/// The generated input files, written once into the target directory.
struct Inputs {
    fastq: PathBuf,
    fasta: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fastq: dir.join("relab-bench.fq"),
        fasta: dir.join("relab-bench.fa"),
    };
    let fastq = generate();
    std::fs::write(&inputs.fasta, to_fasta(&fastq)).unwrap();
    std::fs::write(&inputs.fastq, fastq).unwrap();
    inputs
}

/// A benchmark case: relab's arguments and the equivalent vsearch ones.
struct Case<'a> {
    name: &'static str,
    relab: Vec<&'a str>,
    vsearch: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let case = |name, file, relab: &[&'static str], vsearch: &[&'static str]| {
        let output: &[&str] = if file == fa {
            &["--fastaout", "/dev/null", "--fasta_width", "0"]
        } else {
            &["--fastqout", "/dev/null"]
        };
        Case {
            name,
            relab: [&["-i", file, "-o", "/dev/null"], relab].concat(),
            vsearch: [&["--fastx_filter", file, "--quiet"], output, vsearch].concat(),
        }
    };
    vec![
        case("copy", fq, &[], &[]),
        case(
            "prefix",
            fq,
            &["--relabel", "R", "--sizeout"],
            &["--relabel", "R", "--sizeout"],
        ),
        case("self", fq, &["--relabel-self"], &["--relabel_self"]),
        case("sha1", fq, &["--relabel-sha1"], &["--relabel_sha1"]),
        case(
            "md5_keep",
            fq,
            &["--relabel-md5", "--relabel-keep"],
            &["--relabel_md5", "--relabel_keep"],
        ),
        case(
            "strip",
            fq,
            &["--xsize", "--xee", "--xlength"],
            &["--xsize", "--xee", "--xlength"],
        ),
        case(
            "annotate",
            fq,
            &[
                "--sizeout",
                "--lengthout",
                "--sample",
                "S1",
                "--label-suffix",
                ";",
            ],
            &[
                "--sizeout",
                "--lengthout",
                "--sample",
                "S1",
                "--label_suffix",
                ";",
            ],
        ),
        case("ee", fq, &["--ee"], &["--eeout"]),
        case(
            "fasta_prefix",
            fa,
            &["--relabel", "R", "--sizeout"],
            &["--relabel", "R", "--sizeout"],
        ),
        case(
            "fasta_strip",
            fa,
            &["--xsize", "--xee", "--xlength"],
            &["--xsize", "--xee", "--xlength"],
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
    let relab = env!("CARGO_BIN_EXE_relab");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("relab/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("relab/{}", case.name), |b| {
            b.iter(|| run(relab, &case.relab))
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

criterion_group!(benches, bench_format, bench_binary);
criterion_main!(benches);
