//! Benchmarks of aligning all pairs.
//!
//! `binary` runs the allpr program end to end on 600 sequences of 250
//! bases, copies of 20 templates with one base in twenty changed, inserted
//! or deleted: 179,700 alignments, accepting all, and at --id 0.9. Set
//! `VSEARCH` to a vsearch executable to time `vsearch --allpairs_global
//! --threads 1` on the same file, after checking that both write the same
//! hits:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench allpr
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const SEQUENCES: usize = 600;
const LENGTH: usize = 250;
const TEMPLATES: usize = 20;

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

fn input_file() -> PathBuf {
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
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("allpr-bench.fa");
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

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fa = input.to_str().unwrap();
    let allpr = env!("CARGO_BIN_EXE_allpr");
    let vsearch = std::env::var("VSEARCH").ok();
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group("allpr/binary");
    group.throughput(Throughput::Elements(
        (SEQUENCES * (SEQUENCES - 1) / 2) as u64,
    ));
    group.sample_size(10);
    for (name, rargs, vargs) in [
        ("all", vec!["--accept-all"], vec!["--acceptall"]),
        ("id90", vec!["--id", "0.9"], vec!["--id", "0.9"]),
    ] {
        let vbase = [
            "--allpairs_global",
            fa,
            "--quiet",
            "--notrunclabels",
            "--qmask",
            "none",
            "--threads",
            "1",
        ];
        // both write the same hits
        if let Some(vsearch) = &vsearch {
            let (r, v) = (dir.join("allpr-bench-r.uc"), dir.join("allpr-bench-v.uc"));
            run(
                allpr,
                &[&["-i", fa, "-o", r.to_str().unwrap()][..], &rargs].concat(),
            );
            run(
                vsearch,
                &[&vbase[..], &vargs, &["--uc", v.to_str().unwrap()]].concat(),
            );
            assert!(
                std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
                "allpr and vsearch wrote different hits"
            );
        }
        let rargs = [&["-i", fa, "-o", "/dev/null"][..], &rargs].concat();
        group.bench_function(format!("allpr/{name}"), |b| b.iter(|| run(allpr, &rargs)));
        if let Some(vsearch) = &vsearch {
            let vargs = [&vbase[..], &vargs, &["--uc", "/dev/null"]].concat();
            translated::bench(&mut group, format!("vsearch/{name}"), &vargs);
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &vargs))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
