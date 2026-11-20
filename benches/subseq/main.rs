//! Benchmarks of keeping one region of every read.
//!
//! `binary` runs the subseq program end to end on files of FASTQ and FASTA
//! reads: a region in the middle, and from a start to the end. Set
//! `VSEARCH` to a vsearch executable to time `vsearch --fastx_getsubseq`,
//! selecting every record, on the same files for comparison:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench subseq
//! ```

use std::path::{Path, PathBuf};
use std::process::Command;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

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

/// Generates reads as FASTQ and as FASTA.
fn generate() -> (String, String) {
    let mut rng = Rng(0x5eed);
    let (mut fastq, mut fasta) = (String::new(), String::new());
    for i in 0..READS {
        let seq: String = (0..LENGTH)
            .map(|_| char::from(b"ACGT"[rng.below(4)]))
            .collect();
        let qual: String = (0..LENGTH)
            .map(|_| char::from(b'#' + rng.below(39) as u8))
            .collect();
        fastq.push_str(&format!("@read{i}\n{seq}\n+\n{qual}\n"));
        fasta.push_str(&format!(">read{i}\n{seq}\n"));
    }
    (fastq, fasta)
}

/// The generated input files, written once into the target directory.
struct Inputs {
    fastq: PathBuf,
    fasta: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (fastq, fasta) = generate();
    let inputs = Inputs {
        fastq: dir.join("subseq-bench.fq"),
        fasta: dir.join("subseq-bench.fa"),
    };
    std::fs::write(&inputs.fastq, fastq).unwrap();
    std::fs::write(&inputs.fasta, fasta).unwrap();
    inputs
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let subseq = env!("CARGO_BIN_EXE_subseq");
    let vsearch = std::env::var("VSEARCH").ok();
    let fq = inputs.fastq.to_str().unwrap();
    let fa = inputs.fasta.to_str().unwrap();

    let mut group = c.benchmark_group("subseq/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    let regions: [(&str, &[&str], &[&str]); 2] = [
        (
            "middle",
            &["--start", "50", "--end", "200"],
            &["--subseq_start", "50", "--subseq_end", "200"],
        ),
        ("to_end", &["--start", "20"], &["--subseq_start", "20"]),
    ];
    for (format, file, output) in [
        ("fastq", fq, &["--fastqout", "/dev/null"][..]),
        (
            "fasta",
            fa,
            &["--fastaout", "/dev/null", "--fasta_width", "0"],
        ),
    ] {
        for (region, args, vargs) in regions {
            let name = format!("{region}_{format}");
            let args = [&["-i", file, "-o", "/dev/null"][..], args].concat();
            group.bench_function(format!("subseq/{name}"), |b| b.iter(|| run(subseq, &args)));
            if let Some(vsearch) = &vsearch {
                let args = [
                    &[
                        "--fastx_getsubseq",
                        file,
                        "--quiet",
                        "--notrunclabels",
                        "--label_substr_match",
                        "--label",
                        "",
                    ][..],
                    output,
                    vargs,
                ]
                .concat();
                translated::bench(&mut group, format!("vsearch/{name}"), &args);
                group.bench_function(format!("vsearch/{name}"), |b| {
                    b.iter(|| run(vsearch, &args))
                });
            }
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
