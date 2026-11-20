//! Benchmarks of dereplication.
//!
//! `cluster` measures the in-memory work (`Clusters::add_fastx` and
//! sorting) on parsed records; `binary` runs the derep program end to end
//! on files: FASTA and FASTQ, with and without `--strand`, without sizes,
//! in input order, Phred+64, the best quality, UC and tab-separated tables,
//! identical labels and prefixes. Set `VSEARCH` to a vsearch executable to
//! time `vsearch --fastx_uniques`, `--derep_smallmem`, `--derep_id` and
//! `--derep_prefix` on the same files for comparison, with
//! `--minseqlength 0 --notrunclabels` as derep keeps every length and whole
//! headers:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench derep
//! ```

use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

use rsearch::cluster::{Cluster, Clusters, reverse_complement, sort_by_size};
use rsearch::fastx::{Reader, Record};

/// Reads drawn from `UNIQUES` distinct sequences of `LENGTH` nucleotides.
const READS: usize = 100_000;
const UNIQUES: usize = 2_000;
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

/// Generates reads with a skewed abundance distribution: read i comes from
/// a unique sequence chosen with a bias towards low indices. A tenth of the
/// reads are lower case, and a tenth are reverse complemented.
fn generate(fastq: bool) -> String {
    let mut rng = Rng(0x5eed);
    let uniques: Vec<Vec<u8>> = (0..UNIQUES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut out = String::new();
    for i in 0..READS {
        // the minimum of two uniform draws favours abundant sequences
        let u = rng.below(UNIQUES).min(rng.below(UNIQUES));
        let mut seq = uniques[u].clone();
        match rng.below(10) {
            0 => seq.make_ascii_lowercase(),
            1 => seq = reverse_complement(&seq),
            _ => {}
        }
        let seq = String::from_utf8(seq).unwrap();
        if fastq {
            let qual: String = (0..LENGTH)
                .map(|_| char::from(b'#' + rng.below(39) as u8))
                .collect();
            out.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        } else {
            out.push_str(&format!(">read{i}\n{seq}\n"));
        }
    }
    out
}

fn parse(text: &str) -> Vec<Record> {
    Reader::new(text.as_bytes())
        .records()
        .map(Result::unwrap)
        .collect()
}

fn dereplicate(reads: &[Record], strand: bool) -> Vec<Cluster> {
    let mut clusters = Clusters::default();
    for read in reads {
        clusters.add_fastx(read, strand);
    }
    let mut clusters = clusters.into_clusters();
    sort_by_size(&mut clusters);
    clusters
}

fn bench_cluster(c: &mut Criterion) {
    let fasta_reads = parse(&generate(false));
    let fastq_reads = parse(&generate(true));

    let mut group = c.benchmark_group("derep/cluster");
    group.throughput(Throughput::Elements(READS as u64));
    for (name, reads, strand) in [
        ("fasta", &fasta_reads, false),
        ("fasta_strand", &fasta_reads, true),
        ("fastq", &fastq_reads, false),
        ("fastq_strand", &fastq_reads, true),
    ] {
        group.bench_function(name, |b| b.iter(|| black_box(dereplicate(reads, strand))));
    }
    group.finish();
}

/// Re-encodes Phred+33 FASTQ text as Phred+64.
fn to_phred64(fastq: &str) -> String {
    fastq
        .lines()
        .enumerate()
        .map(|(i, line)| {
            if i % 4 == 3 {
                line.bytes().map(|q| char::from(q + 31)).collect::<String>() + "\n"
            } else {
                format!("{line}\n")
            }
        })
        .collect()
}

/// The generated input files: FASTA, FASTQ, and the FASTQ as Phred+64.
struct Inputs {
    fasta: PathBuf,
    fastq: PathBuf,
    phred64: PathBuf,
}

/// Writes the generated input files once, into the target directory.
fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        fasta: dir.join("derep-bench.fa"),
        fastq: dir.join("derep-bench.fq"),
        phred64: dir.join("derep-bench-64.fq"),
    };
    let fastq = generate(true);
    std::fs::write(&inputs.fasta, generate(false)).unwrap();
    std::fs::write(&inputs.phred64, to_phred64(&fastq)).unwrap();
    std::fs::write(&inputs.fastq, fastq).unwrap();
    inputs
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

/// A benchmark case: derep's arguments and the equivalent vsearch ones.
struct Case<'a> {
    name: &'static str,
    derep: Vec<&'a str>,
    vsearch: Vec<&'a str>,
}

fn cases<'a>(inputs: &'a Inputs) -> Vec<Case<'a>> {
    let fa = inputs.fasta.to_str().unwrap();
    let fq = inputs.fastq.to_str().unwrap();
    let fq64 = inputs.phred64.to_str().unwrap();
    let common = ["--minseqlength", "0", "--notrunclabels", "--quiet"];
    // dereplicating `file` with derep's `derep` arguments, and with vsearch
    // `command` written to `out`, with its `vsearch` arguments
    let case = |name, file, derep: &[&'a str], command, out, vsearch: &[&'a str]| Case {
        name,
        derep: [&["-i", file, "-o", "/dev/null"], derep].concat(),
        vsearch: [&[command, file, out, "/dev/null"], &common[..], vsearch].concat(),
    };
    let uniques = |name, file, derep: &[&'a str], out, vsearch: &[&'a str]| {
        case(name, file, derep, "--fastx_uniques", out, vsearch)
    };
    let (fasta, fastq) = ("--fastaout", "--fastqout");
    let sizes = ["--sizein", "--sizeout"];
    let strand = ["--sizein", "--sizeout", "--strand", "both"];
    vec![
        uniques("fasta", fa, &[], fasta, &sizes),
        uniques("fasta_strand", fa, &["--strand"], fasta, &strand),
        uniques("fastq", fq, &[], fastq, &sizes),
        uniques("fastq_strand", fq, &["--strand"], fastq, &strand),
        uniques("unsized", fa, &["--no-sizein", "--no-sizeout"], fasta, &[]),
        case(
            "unsorted",
            fa,
            &["--unsorted"],
            "--derep_smallmem",
            fasta,
            &sizes,
        ),
        uniques(
            "phred64",
            fq64,
            &["--ascii", "64", "--ascii-out", "64"],
            fastq,
            &[
                "--sizein",
                "--sizeout",
                "--fastq_ascii",
                "64",
                "--fastq_asciiout",
                "64",
            ],
        ),
        uniques(
            "best_qual",
            fq,
            &["--best-qual"],
            fastq,
            &["--sizein", "--sizeout", "--fastq_qout_max"],
        ),
        uniques(
            "uc",
            fa,
            &["--strand", "--uc", "/dev/null", "--tabbedout", "/dev/null"],
            fasta,
            &[
                "--sizein",
                "--sizeout",
                "--strand",
                "both",
                "--uc",
                "/dev/null",
                "--tabbedout",
                "/dev/null",
            ],
        ),
        case(
            "same_label",
            fa,
            &["--same-label"],
            "--derep_id",
            "--output",
            &sizes,
        ),
        case(
            "prefix",
            fa,
            &["--prefix"],
            "--derep_prefix",
            "--output",
            &sizes,
        ),
    ]
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let derep = env!("CARGO_BIN_EXE_derep");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("derep/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("derep/{}", case.name), |b| {
            b.iter(|| run(derep, &case.derep))
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

criterion_group!(benches, bench_cluster, bench_binary);
criterion_main!(benches);
