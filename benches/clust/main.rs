//! Benchmarks of clustering.
//!
//! `binary` runs the clust program end to end on 10,000 sequences of 250
//! bases, copies of 500 templates with one base in twenty changed,
//! inserted or deleted, at --id 0.97, and on both strands at --id 0.9 with
//! --maxaccepts 4. Set `VSEARCH` to a vsearch executable to time `vsearch
//! --cluster_smallmem --usersort --threads 1` on the same file, after
//! checking that both write the same centroids and hits:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench clust
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const SEQUENCES: usize = 10_000;
const LENGTH: usize = 250;
const TEMPLATES: usize = 500;

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
        let size = 1 + rng.below(20);
        fasta.push_str(&format!(
            ">s{i};size={size}\n{}\n",
            String::from_utf8(seq).unwrap()
        ));
    }
    let path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("clust-bench.fa");
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

fn bench_binary(c: &mut Criterion) {
    let input = input_file();
    let fa = input.to_str().unwrap();
    let clust = env!("CARGO_BIN_EXE_clust");
    let vsearch = std::env::var("VSEARCH").ok();
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group("clust/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let wide = ["--id", "0.9", "--maxaccepts", "4", "--strand", "both"];
    for (name, args) in [("id97", &["--id", "0.97"][..]), ("both", &wide[..])] {
        let vbase = [
            "--cluster_smallmem",
            fa,
            "--usersort",
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
        // both write the same centroids and hits
        if let Some(vsearch) = &vsearch {
            let path = |name: &str| dir.join(format!("clust-bench-{name}"));
            let (rc, ru, vc, vu) = (path("r.fa"), path("r.uc"), path("v.fa"), path("v.uc"));
            let s = |p: &PathBuf| p.to_str().unwrap().to_string();
            run(
                clust,
                &[&["-i", fa, "-o", &s(&rc), "--hits", &s(&ru)][..], args].concat(),
            );
            run(
                vsearch,
                &[&vbase[..], args, &["--centroids", &s(&vc), "--uc", &s(&vu)]].concat(),
            );
            for (r, v) in [(rc, vc), (ru, vu)] {
                assert!(
                    std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
                    "clust and vsearch wrote different clusters"
                );
            }
        }
        let rargs = [&["-i", fa, "-o", "/dev/null"][..], args].concat();
        group.bench_function(format!("clust/{name}"), |b| b.iter(|| run(clust, &rargs)));
        if let Some(vsearch) = &vsearch {
            let vargs = [&vbase[..], args, &["--centroids", "/dev/null"]].concat();
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &vargs))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
