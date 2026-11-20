//! Benchmarks of quality conversion.
//!
//! `converter` measures `Converter::convert` on parsed reads; `binary` runs
//! the conv program end to end on files: Phred+33 to Phred+64 and back,
//! Solexa scores, limited output scores, FASTA to FASTQ, SFF to FASTQ with
//! and without clipping, FASTQ and FASTA to wrapped or unwrapped FASTA, and
//! UDB (written by mkudb from the FASTA reads) to FASTA. Set `VSEARCH` to a
//! vsearch executable to time `vsearch --fastq_convert`, `--fasta2fastq`,
//! `--sff_convert`, `--udb2fasta` (after checking that it writes conv's
//! FASTA) and, for FASTA output,
//! `--fastx_subsample --sample_pct 100`, which writes every read, on the
//! same files for comparison:
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

/// Pads `bytes` with zeros to a multiple of 8, as SFF sections are.
fn pad8(bytes: &mut Vec<u8>) {
    bytes.resize(bytes.len().div_ceil(8) * 8, 0);
}

/// Builds an SFF file of the reads of `fastq` (Phred+33), with 400 flows
/// per read, clip points 5 bases in from each end and an index at the end.
fn to_sff(fastq: &str) -> Vec<u8> {
    const FLOWS: u16 = 400;
    let header_length = (31 + FLOWS + 4).div_ceil(8) * 8;
    let lines: Vec<&str> = fastq.lines().collect();
    let mut body = Vec::new();
    for read in lines.chunks(4) {
        let (name, bases, quals) = (&read[0][1..], read[1].as_bytes(), read[3].as_bytes());
        let n = bases.len() as u16;
        let mut header = Vec::new();
        header.extend(
            (16 + name.len() as u16)
                .div_ceil(8)
                .wrapping_mul(8)
                .to_be_bytes(),
        );
        header.extend((name.len() as u16).to_be_bytes());
        header.extend(u32::from(n).to_be_bytes());
        for clip in [5, n - 5, 0, 0] {
            header.extend(clip.to_be_bytes());
        }
        header.extend(name.as_bytes());
        pad8(&mut header);
        body.extend(header);
        let mut data = vec![0; 2 * usize::from(FLOWS)];
        data.extend(vec![1; bases.len()]);
        data.extend(bases);
        data.extend(quals.iter().map(|q| q - 33));
        pad8(&mut data);
        body.extend(data);
    }
    let mut index = b".mft1.00".to_vec();
    index.extend([0; 8]);
    let index_offset = u64::from(header_length) + body.len() as u64;
    let mut header = b".sff".to_vec();
    header.extend(1u32.to_be_bytes());
    header.extend(index_offset.to_be_bytes());
    header.extend((index.len() as u32).to_be_bytes());
    header.extend((lines.len() as u32 / 4).to_be_bytes());
    header.extend(header_length.to_be_bytes());
    header.extend(4u16.to_be_bytes());
    header.extend(FLOWS.to_be_bytes());
    header.push(1);
    header.extend(b"TACG".repeat(usize::from(FLOWS) / 4));
    header.extend(b"TCAG");
    pad8(&mut header);
    pad8(&mut index);
    [header, body, index].concat()
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
    sff: PathBuf,
    udb: PathBuf,
}

fn input_files() -> Inputs {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let inputs = Inputs {
        phred33: dir.join("conv-bench.fq"),
        phred64: dir.join("conv-bench-64.fq"),
        solexa: dir.join("conv-bench-solexa.fq"),
        fasta: dir.join("conv-bench.fa"),
        sff: dir.join("conv-bench.sff"),
        udb: dir.join("conv-bench.udb"),
    };
    let phred33 = generate(33, 0);
    std::fs::write(&inputs.sff, to_sff(&phred33)).unwrap();
    std::fs::write(&inputs.fasta, to_fasta(&phred33)).unwrap();
    std::fs::write(&inputs.phred33, phred33).unwrap();
    std::fs::write(&inputs.phred64, generate(64, 0)).unwrap();
    std::fs::write(&inputs.solexa, generate(64, -5)).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_mkudb"))
        .args(["-i", inputs.fasta.to_str().unwrap()])
        .args(["-o", inputs.udb.to_str().unwrap()])
        .stderr(std::process::Stdio::null())
        .status()
        .unwrap();
    assert!(status.success());
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
    let sff = inputs.sff.to_str().unwrap();
    let udb = inputs.udb.to_str().unwrap();
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
        case("fasta", fa, &[], "--fasta2fastq", &[]),
        case("sff", sff, &[], "--sff_convert", &[]),
        case(
            "sff_clip",
            sff,
            &["--clip"],
            "--sff_convert",
            &["--sff_clip"],
        ),
    ]
    .into_iter()
    .chain(
        [
            ("fastq_to_fasta", fq, &[][..], &[][..]),
            (
                "fastq_to_fasta_unwrapped",
                fq,
                &["--fasta-width", "0"],
                &["--fasta_width", "0"],
            ),
            ("fasta_rewrap", fa, &[], &[]),
        ]
        .map(|(name, file, conv, vsearch)| Case {
            name,
            conv: [&["-i", file, "-o", "/dev/null", "--fasta"], conv].concat(),
            vsearch: [
                &[
                    "--fastx_subsample",
                    file,
                    "--sample_pct",
                    "100",
                    "--minseqlength",
                    "0",
                    "--notrunclabels",
                    "--fastaout",
                    "/dev/null",
                    "--quiet",
                ],
                vsearch,
            ]
            .concat(),
        }),
    )
    .chain([Case {
        name: "udb_to_fasta",
        conv: vec!["-i", udb, "-o", "/dev/null", "--fasta"],
        vsearch: vec!["--udb2fasta", udb, "--output", "/dev/null", "--quiet"],
    }])
    .collect()
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program).args(args).status().unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let inputs = input_files();
    let conv = env!("CARGO_BIN_EXE_conv");
    let vsearch = std::env::var("VSEARCH").ok();

    // conv writes the FASTA vsearch --udb2fasta writes
    if let Some(vsearch) = &vsearch {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
        let (r, v) = (dir.join("conv-bench-r.fa"), dir.join("conv-bench-v.fa"));
        let udb = inputs.udb.to_str().unwrap();
        run(conv, &["-i", udb, "--fasta", "-o", r.to_str().unwrap()]);
        run(
            vsearch,
            &[
                "--udb2fasta",
                udb,
                "--quiet",
                "--output",
                v.to_str().unwrap(),
            ],
        );
        assert!(std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap());
    }

    let mut group = c.benchmark_group("conv/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for case in cases(&inputs) {
        group.bench_function(format!("conv/{}", case.name), |b| {
            b.iter(|| run(conv, &case.conv))
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

criterion_group!(benches, bench_converter, bench_binary);
criterion_main!(benches);
