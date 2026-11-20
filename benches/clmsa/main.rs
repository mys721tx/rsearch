//! Benchmarks of aligning clusters.
//!
//! `binary` runs the clmsa program end to end on 20,000 sequences of 250
//! bases. They are copies of 200 templates with one base in twenty changed,
//! inserted or deleted. Their clusters at --id 0.9 are made by clust
//! beforehand. It writes the multiple alignments, the consensus sequences
//! and the profiles. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --cluster_smallmem --usersort --msaout --consout --profile
//! --threads 1` against clust then clmsa three times, after checking that
//! both write the same:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench clmsa
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const SEQUENCES: usize = 20_000;
const LENGTH: usize = 250;
const TEMPLATES: usize = 200;

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

fn input_file(dir: &Path) -> PathBuf {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut fasta = String::new();
    for i in 0..SEQUENCES {
        let mut seq = Vec::new();
        for &base in &templates[rng.below(TEMPLATES)] {
            match rng.below(60) {
                0 => seq.push(b"ACGT"[rng.below(4)]),
                1 => {}
                2 => seq.extend([base, b"ACGT"[rng.below(4)]]),
                _ => seq.push(base),
            }
        }
        fasta.push_str(&format!(">s{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    let path = dir.join("clmsa-bench.fa");
    std::fs::write(&path, fasta).unwrap();
    path
}

fn run(program: &str, args: &[&str]) {
    let status = Command::new(program)
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success(), "{program} {args:?} failed");
}

fn output(program: &str, args: &[String]) -> Vec<u8> {
    let output = Command::new(program)
        .args(args)
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(output.status.success(), "{program} {args:?} failed");
    output.stdout
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let input = input_file(dir);
    let fa = input.to_str().unwrap();
    let (clust, clmsa) = (env!("CARGO_BIN_EXE_clust"), env!("CARGO_BIN_EXE_clmsa"));
    let vsearch = std::env::var("VSEARCH").ok();
    let hits = dir.join("clmsa-bench.uc");
    let hits = hits.to_str().unwrap();
    let cluster = ["-i", fa, "-o", "/dev/null", "--hits", hits, "--id", "0.9"];
    run(clust, &cluster);
    let formats: [&[&str]; 3] = [&[], &["--consensus"], &["--profile"]];
    let align = |format: &[&'static str]| -> Vec<String> {
        [&["--hits", hits, "-i", fa][..], format]
            .concat()
            .into_iter()
            .map(String::from)
            .collect()
    };
    let path = |name: &str| dir.join(format!("clmsa-bench-{name}"));
    let (msa, cons, prof) = (path("v.msa"), path("v.cons"), path("v.prof"));
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();
    let (msa, cons, prof) = (s(&msa), s(&cons), s(&prof));
    let vargs = [
        "--cluster_smallmem",
        fa,
        "--usersort",
        "--id",
        "0.9",
        "--msaout",
        &msa,
        "--consout",
        &cons,
        "--profile",
        &prof,
        "--sizein",
        "--sizeout",
        "--quiet",
        "--notrunclabels",
        "--qmask",
        "soft",
        "--fasta_width",
        "0",
        "--threads",
        "1",
    ];
    // both write the same
    if let Some(vsearch) = &vsearch {
        run(vsearch, &vargs);
        for (format, written) in formats.iter().zip([&msa, &cons, &prof]) {
            assert!(
                output(clmsa, &align(format)) == std::fs::read(written).unwrap(),
                "clmsa and vsearch wrote different alignments"
            );
        }
    }
    let mut group = c.benchmark_group("clmsa/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    for (name, format) in ["msa", "consensus", "profile"].into_iter().zip(formats) {
        let mut args = align(format);
        args.extend(["-o".to_string(), "/dev/null".to_string()]);
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        group.bench_function(format!("clmsa/{name}"), |b| b.iter(|| run(clmsa, &args)));
    }
    if let Some(vsearch) = &vsearch {
        group.bench_function("clust+clmsa/all", |b| {
            b.iter(|| {
                run(clust, &cluster);
                for format in formats {
                    let mut args = align(format);
                    args.extend(["-o".to_string(), "/dev/null".to_string()]);
                    let args: Vec<&str> = args.iter().map(String::as_str).collect();
                    run(clmsa, &args);
                }
            })
        });
        translated::bench(&mut group, "vsearch/all", &vargs);
        group.bench_function("vsearch/all", |b| b.iter(|| run(vsearch, &vargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
