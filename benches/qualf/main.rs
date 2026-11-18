//! Benchmarks of quality filtering.
//!
//! `filter` measures `Filter::apply` and `label` measures `Label::format`
//! on parsed reads; `binary` runs the qualf program end to end on files:
//! FASTQ filtering, relabelling, paired reads and FASTA input. Set
//! `VSEARCH` to a vsearch executable to time `vsearch --fastq_filter` on
//! the same files for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench qualf
//! ```

use std::hint::black_box;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastx::{Reader, Record};
use rsearch::filter::{Filter, expected_errors};
use rsearch::label::{Label, Relabel, abundance};

/// Reads of `LENGTH` nucleotides.
const READS: usize = 100_000;
const LENGTH: usize = 250;

const FORWARD_SEED: u64 = 0x5eed;
const REVERSE_SEED: u64 = 0xbeef;

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

/// Generates Illumina-like reads: quality starts high and drifts down
/// along the read, with occasional N bases. Different seeds give the
/// forward and reverse reads of pairs.
fn generate(seed: u64) -> String {
    let mut rng = Rng(seed);
    let mut out = String::new();
    for i in 0..READS {
        let mut seq = String::with_capacity(LENGTH);
        let mut qual = String::with_capacity(LENGTH);
        let mut q = 30 + rng.below(12) as i64;
        for _ in 0..LENGTH {
            q = (q + rng.below(4) as i64 - 2).clamp(2, 41);
            seq.push(if rng.below(200) == 0 {
                'N'
            } else {
                b"ACGT"[rng.below(4)] as char
            });
            qual.push((q as u8 + 33) as char);
        }
        out.push_str(&format!(
            "@read{i};size={}\n{seq}\n+\n{qual}\n",
            1 + rng.below(5)
        ));
    }
    out
}

fn parse(text: &str) -> Vec<Record> {
    Reader::new(BufReader::new(text.as_bytes()))
        .records()
        .map(Result::unwrap)
        .collect()
}

/// Typical settings of each kind.
fn filters() -> [(&'static str, Filter); 3] {
    [
        (
            "maxee",
            Filter {
                max_ee: Some(1.0),
                ..Filter::default()
            },
        ),
        (
            "truncate",
            Filter {
                strip_left: 10,
                trunc_qual: Some(20),
                min_len: 100,
                max_ns: Some(0),
                ..Filter::default()
            },
        ),
        (
            "truncee",
            Filter {
                trunc_ee: Some(1.0),
                trunc_len: Some(200),
                min_qual: 5,
                ..Filter::default()
            },
        ),
    ]
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

fn bench_filter(c: &mut Criterion) {
    let reads = parse(&generate(FORWARD_SEED));

    let mut group = c.benchmark_group("qualf/filter");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, filter) in filters() {
        group.bench_function(name, |b| {
            b.iter(|| {
                reads
                    .iter()
                    .filter(|r| {
                        !filter
                            .apply(r.seq(), Some(r.qual()), abundance(r.id()))
                            .unwrap()
                            .discarded
                    })
                    .count()
            })
        });
    }
    group.bench_function("expected_errors", |b| {
        b.iter(|| {
            reads
                .iter()
                .map(|r| black_box(expected_errors(r.qual())))
                .sum::<f64>()
        })
    });
    group.finish();
}

/// Typical label settings of each kind.
fn labels() -> [(&'static str, Label); 3] {
    [
        (
            "number",
            Label {
                relabel: Relabel::Prefix("R".into()),
                sizeout: true,
                eeout: true,
                ..Label::default()
            },
        ),
        (
            "sha1",
            Label {
                relabel: Relabel::Sha1,
                lengthout: true,
                ..Label::default()
            },
        ),
        (
            "md5",
            Label {
                relabel: Relabel::Md5,
                ..Label::default()
            },
        ),
    ]
}

fn bench_label(c: &mut Criterion) {
    let reads = parse(&generate(FORWARD_SEED));

    let mut group = c.benchmark_group("qualf/label");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, label) in labels() {
        group.bench_function(name, |b| {
            b.iter(|| {
                for (i, r) in reads.iter().enumerate() {
                    let ee = Some(expected_errors(r.qual()));
                    black_box(label.format(r.header(), r.seq(), 1, i as u64 + 1, ee));
                }
            })
        });
    }
    group.finish();
}

/// The generated input files, written once into the target directory.
struct Inputs {
    forward: PathBuf,
    reverse: PathBuf,
    fasta: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        forward: dir.join("qualf-bench.fq"),
        reverse: dir.join("qualf-bench-rev.fq"),
        fasta: dir.join("qualf-bench.fa"),
    };
    let forward = generate(FORWARD_SEED);
    std::fs::write(&inputs.fasta, to_fasta(&forward)).unwrap();
    std::fs::write(&inputs.forward, forward).unwrap();
    std::fs::write(&inputs.reverse, generate(REVERSE_SEED)).unwrap();
    inputs
}

/// A benchmark case: qualf's arguments and the equivalent vsearch ones,
/// each after the input and output files.
struct Case {
    name: &'static str,
    qualf: Vec<String>,
    vsearch: Vec<String>,
}

fn args(list: &[&str]) -> Vec<String> {
    list.iter().map(|s| s.to_string()).collect()
}

fn cases(inputs: &Inputs) -> Vec<Case> {
    let fq = inputs.forward.to_str().unwrap();
    let rev = inputs.reverse.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();
    // filtering the forward reads, written as FASTQ
    let fastq = |name, qualf: &[&str], vsearch: &[&str]| Case {
        name,
        qualf: args(&[&["-i", fq, "-o", "/dev/null"], qualf].concat()),
        vsearch: args(
            &[
                &["--fastq_filter", fq, "--fastqout", "/dev/null", "--quiet"],
                vsearch,
            ]
            .concat(),
        ),
    };
    vec![
        fastq("maxee", &["--prob", "1"], &["--fastq_maxee", "1"]),
        fastq(
            "truncate",
            &[
                "--strip-left",
                "10",
                "--trunc-qual",
                "20",
                "--min",
                "100",
                "--bad",
                "0",
            ],
            &[
                "--fastq_stripleft",
                "10",
                "--fastq_truncqual",
                "20",
                "--fastq_minlen",
                "100",
                "--fastq_maxns",
                "0",
            ],
        ),
        fastq(
            "truncee",
            &["--trunc-prob", "1", "--trunc-len", "200", "--min-qual", "5"],
            &[
                "--fastq_truncee",
                "1",
                "--fastq_trunclen",
                "200",
                "--fastq_minqual",
                "5",
            ],
        ),
        fastq(
            "relabel",
            &["--prob", "1", "--relabel", "R", "--sizeout", "--ee"],
            &[
                "--fastq_maxee",
                "1",
                "--relabel",
                "R",
                "--sizeout",
                "--eeout",
            ],
        ),
        fastq(
            "sha1",
            &["--relabel-sha1", "--lengthout"],
            &["--relabel_sha1", "--lengthout"],
        ),
        fastq("md5", &["--relabel-md5"], &["--relabel_md5"]),
        Case {
            name: "paired",
            qualf: args(&[
                "-i",
                fq,
                "--reverse",
                rev,
                "-o",
                "/dev/null",
                "--rev-out",
                "/dev/null",
                "--prob",
                "1",
            ]),
            vsearch: args(&[
                "--fastq_filter",
                fq,
                "--reverse",
                rev,
                "--fastqout",
                "/dev/null",
                "--fastqout_rev",
                "/dev/null",
                "--quiet",
                "--fastq_maxee",
                "1",
            ]),
        },
        Case {
            name: "fasta",
            qualf: args(&[
                "-i",
                fa,
                "-o",
                "/dev/null",
                "--min",
                "100",
                "--bad",
                "0",
                "--min-size",
                "2",
            ]),
            vsearch: args(&[
                "--fastx_filter",
                fa,
                "--fastaout",
                "/dev/null",
                "--quiet",
                "--fasta_width",
                "0",
                "--fastq_minlen",
                "100",
                "--fastq_maxns",
                "0",
                "--minsize",
                "2",
            ]),
        },
    ]
}

fn run(program: &str, args: &[String]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let qualf = env!("CARGO_BIN_EXE_qualf");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("qualf/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("qualf/{}", case.name), |b| {
            b.iter(|| run(qualf, &case.qualf))
        });
        if let Some(vsearch) = &vsearch {
            group.bench_function(format!("vsearch/{}", case.name), |b| {
                b.iter(|| run(vsearch, &case.vsearch))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_filter, bench_label, bench_binary);
criterion_main!(benches);
