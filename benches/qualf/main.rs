//! Benchmarks of quality filtering.
//!
//! `filter` measures `Filter::apply` on parsed reads; `binary` runs the
//! qualf program end to end on a file. Set `VSEARCH` to a vsearch
//! executable to time `vsearch --fastq_filter` on the same file for
//! comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench qualf
//! ```

use std::hint::black_box;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::Command;

use bio::io::fastq::Record;
use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::fastq::Reader;
use rsearch::filter::{Filter, abundance, expected_errors};

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

/// Generates Illumina-like reads: quality starts high and drifts down
/// along the read, with occasional N bases.
fn generate() -> String {
    let mut rng = Rng(0x5eed);
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

fn bench_filter(c: &mut Criterion) {
    let reads = parse(&generate());

    let mut group = c.benchmark_group("qualf/filter");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, filter) in filters() {
        group.bench_function(name, |b| {
            b.iter(|| {
                reads
                    .iter()
                    .filter(|r| !filter.apply(r.seq(), r.qual(), abundance(r.id())).discarded)
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

/// Writes the generated input file once, into the target directory.
fn input_file() -> PathBuf {
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("qualf-bench.fq");
    std::fs::write(&path, generate()).unwrap();
    path
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn bench_binary(c: &mut Criterion) {
    let file = input_file();
    let file = file.to_str().unwrap();
    let qualf = env!("CARGO_BIN_EXE_qualf");
    let vsearch = std::env::var("VSEARCH").ok();

    let cases: [(&str, &[&str], &[&str]); 3] = [
        ("maxee", &["--prob", "1"], &["--fastq_maxee", "1"]),
        (
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
        (
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
    ];

    let mut group = c.benchmark_group("qualf/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, qualf_args, vsearch_args) in cases {
        let mut args = vec!["-i", file, "-o", "/dev/null"];
        args.extend(qualf_args);
        group.bench_function(format!("qualf/{name}"), |b| b.iter(|| run(qualf, &args)));
        if let Some(vsearch) = &vsearch {
            let mut args = vec!["--fastq_filter", file, "--fastqout", "/dev/null", "--quiet"];
            args.extend(vsearch_args);
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_filter, bench_binary);
criterion_main!(benches);
