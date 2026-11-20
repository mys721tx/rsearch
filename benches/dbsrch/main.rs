//! Benchmarks of searching a database.
//!
//! `binary` runs the dbsrch program end to end: 2,000 queries of 250
//! bases against a database of 10,000, copies of 200 templates with one
//! base in twenty changed, inserted or deleted, at the defaults with --id
//! 0.97, and on both strands with --id 0.9, --maxaccepts 8 and --maxrejects
//! 64. Set `VSEARCH` to a vsearch executable to time `vsearch
//! --usearch_global --threads 1` on the same files, after checking that
//! both write the same hits:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench dbsrch
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const QUERIES: usize = 2_000;
const TARGETS: usize = 10_000;
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

/// the queries and the database.
fn input_files() -> (PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut fasta = |prefix: &str, n: usize| {
        let mut text = String::new();
        for i in 0..n {
            let mut seq = Vec::new();
            for &base in &templates[rng.below(TEMPLATES)] {
                match rng.below(60) {
                    0 => seq.push(b"ACGT"[rng.below(4)]),
                    1 => {}
                    2 => seq.extend([base, b"ACGT"[rng.below(4)]]),
                    _ => seq.push(base),
                }
            }
            text.push_str(&format!(
                ">{prefix}{i}\n{}\n",
                String::from_utf8(seq).unwrap()
            ));
        }
        text
    };
    let db = fasta("t", TARGETS);
    let q = fasta("q", QUERIES);
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let (qp, dbp) = (
        dir.join("dbsrch-bench-q.fa"),
        dir.join("dbsrch-bench-db.fa"),
    );
    std::fs::write(&qp, q).unwrap();
    std::fs::write(&dbp, db).unwrap();
    (qp, dbp)
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
    let (q, db) = input_files();
    let (q, db) = (q.to_str().unwrap(), db.to_str().unwrap());
    let dbsrch = env!("CARGO_BIN_EXE_dbsrch");
    let vsearch = std::env::var("VSEARCH").ok();
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"));
    let mut group = c.benchmark_group("dbsrch/binary");
    group.throughput(Throughput::Elements(QUERIES as u64));
    group.sample_size(10);
    let wide = [
        "--id",
        "0.9",
        "--maxaccepts",
        "8",
        "--maxrejects",
        "64",
        "--strand",
        "both",
    ];
    for (name, args) in [("id97", &["--id", "0.97"][..]), ("both", &wide[..])] {
        let vbase = [
            "--usearch_global",
            q,
            "--db",
            db,
            "--quiet",
            "--notrunclabels",
            "--qmask",
            "soft",
            "--dbmask",
            "soft",
            "--uc_allhits",
            "--threads",
            "1",
        ];
        // both write the same hits
        if let Some(vsearch) = &vsearch {
            let (r, v) = (dir.join("dbsrch-bench-r.uc"), dir.join("dbsrch-bench-v.uc"));
            run(
                dbsrch,
                &[&["-i", q, "--db", db, "-o", r.to_str().unwrap()][..], args].concat(),
            );
            run(
                vsearch,
                &[&vbase[..], args, &["--uc", v.to_str().unwrap()]].concat(),
            );
            assert!(
                std::fs::read(&r).unwrap() == std::fs::read(&v).unwrap(),
                "dbsrch and vsearch wrote different hits"
            );
        }
        let rargs = [&["-i", q, "--db", db, "-o", "/dev/null"][..], args].concat();
        group.bench_function(format!("dbsrch/{name}"), |b| b.iter(|| run(dbsrch, &rargs)));
        if let Some(vsearch) = &vsearch {
            let vargs = [&vbase[..], args, &["--uc", "/dev/null"]].concat();
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
