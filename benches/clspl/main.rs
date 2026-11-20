//! Benchmarks of splitting clusters.
//!
//! `binary` runs the clspl program end to end on 20,000 sequences of 250
//! bases. They are copies of 200 templates with one base in twenty changed,
//! inserted or deleted. Their clusters at --id 0.97 are made by clust
//! beforehand. Set `VSEARCH` to a vsearch executable to time `vsearch
//! --cluster_smallmem --usersort --clusters --threads 1` against clust then
//! clspl, after checking that both write the same files:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench clspl
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
    let path = dir.join("clspl-bench.fa");
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

/// the files written with `prefix`, in cluster order.
fn files(prefix: &Path) -> Vec<Vec<u8>> {
    (0..)
        .map_while(|i| std::fs::read(format!("{}{i}", prefix.display())).ok())
        .collect()
}

#[path = "../common/translated.rs"]
mod translated;

fn bench_binary(c: &mut Criterion) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let input = input_file(dir);
    let fa = input.to_str().unwrap();
    let (clust, clspl) = (env!("CARGO_BIN_EXE_clust"), env!("CARGO_BIN_EXE_clspl"));
    let vsearch = std::env::var("VSEARCH").ok();
    let hits = dir.join("clspl-bench.uc");
    let hits = hits.to_str().unwrap();
    run(
        clust,
        &["-i", fa, "-o", "/dev/null", "--hits", hits, "--id", "0.97"],
    );
    let (rdir, vdir) = (dir.join("clspl-bench-r"), dir.join("clspl-bench-v"));
    for d in [&rdir, &vdir] {
        let _ = std::fs::remove_dir_all(d);
        std::fs::create_dir_all(d).unwrap();
    }
    let (rprefix, vprefix) = (rdir.join("c"), vdir.join("c"));
    let (rp, vp) = (rprefix.to_str().unwrap(), vprefix.to_str().unwrap());
    let split = ["--hits", hits, "-i", fa, "--prefix", rp];
    let vargs = [
        "--cluster_smallmem",
        fa,
        "--usersort",
        "--id",
        "0.97",
        "--clusters",
        vp,
        "--quiet",
        "--notrunclabels",
        "--qmask",
        "soft",
        "--fasta_width",
        "0",
        "--threads",
        "1",
    ];
    // both write the same files
    if let Some(vsearch) = &vsearch {
        run(clspl, &split);
        run(vsearch, &vargs);
        assert!(
            files(&rprefix) == files(&vprefix),
            "clspl and vsearch wrote different files"
        );
    }
    let mut group = c.benchmark_group("clspl/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    group.bench_function("clspl/split", |b| b.iter(|| run(clspl, &split)));
    if let Some(vsearch) = &vsearch {
        let cluster = ["-i", fa, "-o", "/dev/null", "--hits", hits, "--id", "0.97"];
        group.bench_function("clust+clspl/cluster", |b| {
            b.iter(|| {
                run(clust, &cluster);
                run(clspl, &split);
            })
        });
        translated::bench(&mut group, "vsearch/cluster", &vargs);
        group.bench_function("vsearch/cluster", |b| b.iter(|| run(vsearch, &vargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
