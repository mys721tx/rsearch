//! Benchmarks of SINTAX classification.
//!
//! `binary` runs the sintax program end to end: 10,000 reads of 250 bases
//! with one error in fifty, cut from 50 templates of 1,400 bases, half of
//! them reverse-complemented, against a database of 5,000 mutated copies of
//! the templates, each annotated with its template's taxonomy, at the
//! default word length of 8, searching both strands, from a FASTA file and
//! from mkudb's UDB file of it. Set `VSEARCH` to a vsearch executable to
//! time `vsearch --sintax --threads 1` on the same files, after checking
//! that both write the same table:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench sintax
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};

const READS: usize = 10_000;
const READ_LENGTH: usize = 250;
const DATABASE: usize = 5_000;
const LENGTH: usize = 1_400;
const TEMPLATES: usize = 50;

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

fn mutate(rng: &mut Rng, seq: &mut [u8], one_in: usize) {
    for base in seq.iter_mut() {
        if rng.below(one_in) == 0 {
            *base = b"ACGT"[rng.below(4)];
        }
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn tmp(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

/// writes the reads, the FASTA database and its UDB file.
fn input_files() -> (PathBuf, PathBuf, PathBuf) {
    let mut rng = Rng(0x5eed);
    let templates: Vec<Vec<u8>> = (0..TEMPLATES)
        .map(|_| (0..LENGTH).map(|_| b"ACGT"[rng.below(4)]).collect())
        .collect();
    let mut db = String::new();
    for i in 0..DATABASE {
        let t = rng.below(TEMPLATES);
        let mut seq = templates[t].clone();
        mutate(&mut rng, &mut seq, 20);
        db.push_str(&format!(
            ">d{i};tax=d:D{},p:P{},g:G{t}\n{}\n",
            t % 2,
            t % 7,
            String::from_utf8(seq).unwrap()
        ));
    }
    let mut reads = String::new();
    for i in 0..READS {
        let template = &templates[rng.below(TEMPLATES)];
        let start = rng.below(LENGTH - READ_LENGTH);
        let mut seq = template[start..start + READ_LENGTH].to_vec();
        mutate(&mut rng, &mut seq, 50);
        if rng.below(2) == 0 {
            seq = reverse_complement(&seq);
        }
        reads.push_str(&format!(">r{i}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    let (fa, dbfa, udb) = (
        tmp("sintax-bench.fa"),
        tmp("sintax-bench-db.fa"),
        tmp("sintax-bench-db.udb"),
    );
    std::fs::write(&fa, reads).unwrap();
    std::fs::write(&dbfa, db).unwrap();
    run(
        env!("CARGO_BIN_EXE_mkudb"),
        &["-i", dbfa.to_str().unwrap(), "-o", udb.to_str().unwrap()],
    );
    (fa, dbfa, udb)
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
    let (fa, dbfa, udb) = input_files();
    let fa = fa.to_str().unwrap();
    let sintax = env!("CARGO_BIN_EXE_sintax");
    let vsearch = std::env::var("VSEARCH").ok();

    let mut group = c.benchmark_group("sintax/binary");
    group.throughput(Throughput::Elements(READS as u64));
    group.sample_size(10);
    for (name, db) in [("fasta", &dbfa), ("udb", &udb)] {
        let db = db.to_str().unwrap();
        let (r, v) = (tmp("sintax-bench-r.tsv"), tmp("sintax-bench-v.tsv"));
        let (r, v) = (r.to_str().unwrap(), v.to_str().unwrap());
        let rargs = [
            "-i", fa, "--db", db, "-o", r, "--seed", "1", "--strand", "both",
        ];
        let vargs = [
            "--sintax",
            fa,
            "--db",
            db,
            "--tabbedout",
            v,
            "--randseed",
            "1",
            "--strand",
            "both",
            "--quiet",
            "--notrunclabels",
            "--threads",
            "1",
        ];
        // both write the same table
        if let Some(vsearch) = &vsearch {
            run(sintax, &rargs);
            run(vsearch, &vargs);
            assert!(
                std::fs::read(r).unwrap() == std::fs::read(v).unwrap(),
                "sintax and vsearch wrote different tables"
            );
        }
        group.bench_function(format!("sintax/{name}"), |b| b.iter(|| run(sintax, &rargs)));
        if let Some(vsearch) = &vsearch {
            group.bench_function(format!("vsearch/{name}"), |b| {
                b.iter(|| run(vsearch, &vargs))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
