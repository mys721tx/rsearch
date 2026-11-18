//! Benchmarks of dereplication.
//!
//! `cluster` measures the in-memory work (`Clusters::add_fastx` and
//! sorting) on parsed records; `binary` runs the derep program end to end
//! on files: FASTA and FASTQ, with and without `--strand`, relabelling,
//! whole labels without sizes, and FASTQ written as FASTA.
//! Set `VSEARCH` to a vsearch executable to time `vsearch --fastx_uniques`
//! on the same files for comparison:
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

/// Writes the generated input files once, into the target directory.
fn input_files() -> (PathBuf, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fa, fq) = (dir.join("derep-bench.fa"), dir.join("derep-bench.fq"));
    std::fs::write(&fa, generate(false)).unwrap();
    std::fs::write(&fq, generate(true)).unwrap();
    (fa, fq)
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

/// A benchmark case: the input file, derep's arguments and the equivalent
/// vsearch ones, each after the input and output files.
struct Case<'a> {
    name: &'static str,
    file: &'a str,
    derep: Vec<&'a str>,
    vsearch: Vec<&'a str>,
}

fn cases<'a>(fa: &'a str, fq: &'a str) -> Vec<Case<'a>> {
    let case = |name, file, derep: &[&'a str], vsearch: &[&'a str]| Case {
        name,
        file,
        derep: derep.to_vec(),
        vsearch: vsearch.to_vec(),
    };
    let sizes = ["--sizein", "--sizeout"];
    let strand = ["--sizein", "--sizeout", "--strand", "both"];
    vec![
        case("fasta", fa, &[], &sizes),
        case("fasta_strand", fa, &["--strand"], &strand),
        case("fastq", fq, &[], &sizes),
        case("fastq_strand", fq, &["--strand"], &strand),
        case(
            "relabel",
            fa,
            &["--relabel-sha1", "--lengthout"],
            &["--sizein", "--sizeout", "--relabel_sha1", "--lengthout"],
        ),
        case(
            "unsized",
            fa,
            &["--no-sizein", "--no-sizeout", "--full-labels"],
            &["--notrunclabels"],
        ),
        case("fastq_as_fasta", fq, &["--fasta"], &sizes),
    ]
}

fn bench_binary(c: &mut Criterion) {
    let (fa, fq) = input_files();
    let (fa, fq) = (fa.to_str().unwrap(), fq.to_str().unwrap());
    let derep = env!("CARGO_BIN_EXE_derep");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("derep/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(fa, fq) {
        let mut args = vec!["-i", case.file, "-o", "/dev/null"];
        args.extend(&case.derep);
        group.bench_function(format!("derep/{}", case.name), |b| {
            b.iter(|| run(derep, &args))
        });
        if let Some(vsearch) = &vsearch {
            // FASTQ input is written as FASTQ unless derep is told otherwise
            let out = if case.file == fq && !case.derep.contains(&"--fasta") {
                "--fastqout"
            } else {
                "--fastaout"
            };
            let mut args = vec![
                "--fastx_uniques",
                case.file,
                out,
                "/dev/null",
                "--minseqlength",
                "1",
                "--quiet",
            ];
            args.extend(&case.vsearch);
            group.bench_function(format!("vsearch/{}", case.name), |b| {
                b.iter(|| run(vsearch, &args))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_cluster, bench_binary);
criterion_main!(benches);
