//! Benchmarks of read filtering.
//!
//! `filter` measures `Filter::apply` with filtering settings on parsed
//! reads; `binary` runs the sqfilt program end to end on files: no filter,
//! length, Ns, expected errors and their rate, lowest score, abundance, all
//! of them together, pairs, the share of upper-case residues and FASTA.
//! Set `VSEARCH` to a vsearch executable to time `vsearch --fastx_filter`
//! (or, for the masked share, `--fastx_mask --qmask soft`) with the same
//! filters on the same files for comparison. Both write their discarded
//! reads too.
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench sqfilt
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

/// Generates Illumina-like Phred+33 reads with scores drifting from high to
/// low along the read. `tag` ends each label so that mates can share
/// numbers.
fn generate(tag: &str, seed: u64) -> String {
    let mut rng = Rng(seed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        for _ in 0..LENGTH {
            q = (q + rng.below(4) as i64 - 2).clamp(2, 41);
            seq.push(char::from(b"ACGTACGTACGTACGTN"[rng.below(17)]));
            qual.push(char::from((q + 33) as u8));
        }
        // some reads carry an abundance and a soft-masked stretch
        if i % 3 == 0 {
            seq.replace_range(20..60, &seq[20..60].to_ascii_lowercase());
        }
        let size = 1 + rng.below(5);
        out.push_str(&format!("@read{i}{tag};size={size}\n{seq}\n+\n{qual}\n"));
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

fn bench_filter(c: &mut Criterion) {
    let reads = parse(&generate("", 0x5eed));
    let mut group = c.benchmark_group("sqfilt/filter");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, filter) in [
        (
            "max_ee",
            Filter {
                max_ee: Some(1.0),
                ..Filter::default()
            },
        ),
        (
            "min_qual",
            Filter {
                min_qual: 10,
                ..Filter::default()
            },
        ),
        (
            "all",
            Filter {
                min_len: 100,
                max_ns: Some(2),
                max_ee: Some(2.0),
                max_rate: Some(0.01),
                min_qual: 5,
                ..Filter::default()
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
    reverse: PathBuf,
    fasta: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fastq: dir.join("sqfilt-bench.fq"),
        reverse: dir.join("sqfilt-bench-rev.fq"),
        fasta: dir.join("sqfilt-bench.fa"),
    };
    let fastq = generate("", 0x5eed);
    std::fs::write(&inputs.fasta, to_fasta(&fastq)).unwrap();
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.reverse, generate("/2", 0xfeed)).unwrap();
    inputs
}

/// A benchmark case: sqfilt's arguments and the equivalent vsearch ones.
struct Case<'a> {
    name: &'static str,
    sqfilt: Vec<&'a str>,
    vsearch: Vec<&'a str>,
}

fn cases(inputs: &Inputs) -> Vec<Case<'_>> {
    let fq = inputs.fastq.to_str().unwrap();
    let rev = inputs.reverse.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    let null = "/dev/null";
    let case = |name, file, sqfilt: &[&'static str], vsearch: &[&'static str]| {
        let output: &[&str] = if file == fa {
            &[
                "--fastaout",
                null,
                "--fastaout_discarded",
                null,
                "--fasta_width",
                "0",
            ]
        } else {
            &["--fastqout", null, "--fastqout_discarded", null]
        };
        Case {
            name,
            sqfilt: [&["-i", file, "-o", null, "--discarded", null], sqfilt].concat(),
            vsearch: [
                &["--fastx_filter", file, "--quiet", "--notrunclabels"],
                output,
                vsearch,
            ]
            .concat(),
        }
    };
    vec![
        case("none", fq, &[], &[]),
        case(
            "min_len",
            fq,
            &["--min-len", "240"],
            &["--fastq_minlen", "240"],
        ),
        case("max_ns", fq, &["--max-ns", "10"], &["--fastq_maxns", "10"]),
        case("max_ee", fq, &["--max-ee", "1"], &["--fastq_maxee", "1"]),
        case(
            "max_ee_rate",
            fq,
            &["--max-ee-rate", "0.005"],
            &["--fastq_maxee_rate", "0.005"],
        ),
        case(
            "min_qual",
            fq,
            &["--min-qual", "10"],
            &["--fastq_minqual", "10"],
        ),
        case(
            "size",
            fq,
            &["--min-size", "2", "--max-size", "4"],
            &["--minsize", "2", "--maxsize", "4"],
        ),
        case(
            "all",
            fq,
            &[
                "--min-len",
                "100",
                "--max-ns",
                "20",
                "--max-ee",
                "3",
                "--min-qual",
                "5",
                "--min-size",
                "2",
            ],
            &[
                "--fastq_minlen",
                "100",
                "--fastq_maxns",
                "20",
                "--fastq_maxee",
                "3",
                "--fastq_minqual",
                "5",
                "--minsize",
                "2",
            ],
        ),
        Case {
            name: "pairs",
            sqfilt: vec![
                "-i",
                fq,
                "--reverse",
                rev,
                "-o",
                null,
                "--rev-out",
                null,
                "--discarded",
                null,
                "--rev-discarded",
                null,
                "--max-ee",
                "1",
            ],
            vsearch: vec![
                "--fastx_filter",
                fq,
                "--reverse",
                rev,
                "--quiet",
                "--notrunclabels",
                "--fastqout",
                null,
                "--fastqout_rev",
                null,
                "--fastqout_discarded",
                null,
                "--fastqout_discarded_rev",
                null,
                "--fastq_maxee",
                "1",
            ],
        },
        Case {
            name: "upper_pct",
            sqfilt: vec!["-i", fq, "-o", null, "--min-upper-pct", "90"],
            vsearch: vec![
                "--fastx_mask",
                fq,
                "--quiet",
                "--notrunclabels",
                "--qmask",
                "soft",
                "--fastqout",
                null,
                "--min_unmasked_pct",
                "90",
            ],
        },
        case(
            "fasta",
            fa,
            &["--min-len", "100", "--max-ns", "10"],
            &["--fastq_minlen", "100", "--fastq_maxns", "10"],
        ),
    ]
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let sqfilt = env!("CARGO_BIN_EXE_sqfilt");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("sqfilt/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("sqfilt/{}", case.name), |b| {
            b.iter(|| run(sqfilt, &case.sqfilt))
        });
        if let Some(vsearch) = &vsearch {
            group.bench_function(format!("vsearch/{}", case.name), |b| {
                b.iter(|| run(vsearch, &case.vsearch))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_filter, bench_binary);
criterion_main!(benches);
