//! Benchmarks of chimera detection.
//!
//! `binary` runs the uchime program end to end: de novo on 4,000
//! sequences of 250 bases, sorted by decreasing abundance, abundant
//! variants of 100 parents and rarer two-parent chimeras of them, one base
//! in a hundred changed; and the same sequences but the parents against a
//! database of the parents. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --uchime_denovo` and `--uchime_ref` with `--threads 1` on the
//! same files, after checking that both write the same tables:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench uchime
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const SEQUENCES: usize = 4_000;
const LENGTH: usize = 250;
const PARENTS: usize = 100;

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

/// the input, sorted by decreasing abundance, the queries (the input but
/// the parents) and the parents.
fn input_files(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let parents: Vec<Vec<u8>> = (0..PARENTS)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let (mut all, mut queries, mut db) = (String::new(), String::new(), String::new());
    for (i, p) in parents.iter().enumerate() {
        let record = format!(
            ">p{i};size={}\n{}\n",
            100_000 - i,
            String::from_utf8_lossy(p)
        );
        all.push_str(&record);
        db.push_str(&record);
    }
    for i in PARENTS..SEQUENCES {
        let a = rng.below(PARENTS);
        let mut seq = parents[a].clone();
        // a third of them chimeras
        if rng.below(3) == 0 {
            let b = rng.below(PARENTS);
            let cut = 50 + rng.below(LENGTH - 100);
            seq[cut..].copy_from_slice(&parents[b][cut..]);
        }
        for base in seq.iter_mut() {
            if rng.below(100) == 0 {
                *base = b"ACGT"[rng.below(4)];
            }
        }
        // decreasing abundances
        let record = format!(
            ">s{i};size={}\n{}\n",
            50_000 - i,
            String::from_utf8_lossy(&seq)
        );
        all.push_str(&record);
        queries.push_str(&record);
    }
    let paths = ["all", "queries", "db"].map(|n| dir.join(format!("uchime-bench-{n}.fa")));
    for (path, text) in paths.iter().zip([&all, &queries, &db]) {
        std::fs::write(path, text).unwrap();
    }
    let [a, q, d] = paths;
    (a, q, d)
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
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (all, queries, db) = input_files(dir);
    let s = |p: &PathBuf| p.to_str().unwrap().to_string();
    let (all, queries, db) = (s(&all), s(&queries), s(&db));
    let uchime = env!("CARGO_BIN_EXE_uchime");
    let vsearch = std::env::var("VSEARCH").ok();
    let mut group = c.benchmark_group("uchime/binary");
    group.throughput(Throughput::Elements(SEQUENCES as u64));
    group.sample_size(10);
    let common = [
        "--quiet",
        "--notrunclabels",
        "--qmask",
        "soft",
        "--fasta_width",
        "0",
        "--threads",
        "1",
    ];
    let cases: [(&str, Vec<&str>, Vec<&str>); 2] = [
        ("denovo", vec!["-i", &all], vec!["--uchime_denovo", &all]),
        (
            "ref",
            vec!["-i", &queries, "--db", &db],
            vec!["--uchime_ref", &queries, "--db", &db, "--dbmask", "soft"],
        ),
    ];
    for (name, rargs, vargs) in cases {
        // both write the same tables
        if let Some(vsearch) = &vsearch {
            let (r, v) = (
                dir.join("uchime-bench-r.tab"),
                dir.join("uchime-bench-v.tab"),
            );
            run(
                uchime,
                &[&rargs[..], &["-o", "/dev/null", "--table", &s(&r)]].concat(),
            );
            run(
                vsearch,
                &[&vargs[..], &common, &["--uchimeout", &s(&v)]].concat(),
            );
            assert!(
                std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
                "uchime and vsearch wrote different tables"
            );
        }
        let ra = [&rargs[..], &["-o", "/dev/null"]].concat();
        group.bench_function(format!("uchime/{name}"), |b| b.iter(|| run(uchime, &ra)));
        if let Some(vsearch) = &vsearch {
            let va = [&vargs[..], &common, &["--nonchimeras", "/dev/null"]].concat();
            translated::bench(&mut group, format!("vsearch/{name}"), &va);
            group.bench_function(format!("vsearch/{name}"), |b| b.iter(|| run(vsearch, &va)));
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
