//! Benchmarks of reporting a UDB file's word statistics.
//!
//! `binary` runs the udbstat program end to end on the UDB file mkudb
//! writes for 20,000 random sequences of 1,400 bases at word length 12, a
//! file of about 170 MB. Set `VSEARCH` to a vsearch executable to time
//! `vsearch --udbstats --log` on the same file (vsearch analyses nothing
//! without a log), after checking that both report the same statistics:
//!
//! ```sh
//! VSEARCH=/path/to/vsearch cargo bench --bench udbstat
//! ```

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::{Criterion, criterion_group, criterion_main};

const SEQUENCES: usize = 20_000;
const LENGTH: usize = 1_400;

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
}

fn tmp(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(name)
}

/// writes the FASTA input and builds its UDB file with mkudb.
fn input_file() -> PathBuf {
    let mut rng = Rng(0x5eed);
    let mut fasta = String::new();
    for i in 0..SEQUENCES {
        let seq: String = (0..LENGTH)
            .map(|_| b"ACGT"[(rng.next() % 4) as usize] as char)
            .collect();
        fasta.push_str(&format!(">seq{i}\n{seq}\n"));
    }
    let (fa, udb) = (tmp("udbstat-bench.fa"), tmp("udbstat-bench.udb"));
    std::fs::write(&fa, fasta).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_mkudb"))
        .args([
            "--word-len",
            "12",
            "-i",
            fa.to_str().unwrap(),
            "-o",
            udb.to_str().unwrap(),
        ])
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(status.success());
    udb
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
    let udb = input_file();
    let udb = udb.to_str().unwrap();
    let udbstat = env!("CARGO_BIN_EXE_udbstat");
    let vsearch = std::env::var("VSEARCH").ok();

    // both report the same statistics
    let log = tmp("udbstat-bench.log");
    let vargs = [
        "--udbstats",
        udb,
        "--quiet",
        "--threads",
        "1",
        "--log",
        log.to_str().unwrap(),
    ];
    if let Some(vsearch) = &vsearch {
        let ours = Command::new(udbstat).args(["-i", udb]).output().unwrap();
        let _ = std::fs::remove_file(&log);
        run(vsearch, &vargs);
        let log = std::fs::read_to_string(&log).unwrap();
        // the report: past the log's three opening lines, up to its blank
        // line and three closing ones
        let lines: Vec<&str> = log.lines().collect();
        let theirs: String = lines[3..lines.len() - 4]
            .iter()
            .map(|l| format!("{l}\n"))
            .collect();
        assert_eq!(String::from_utf8(ours.stdout).unwrap(), theirs);
    }

    let mut group = c.benchmark_group("udbstat/binary");
    group.sample_size(10);
    group.bench_function("udbstat", |b| b.iter(|| run(udbstat, &["-i", udb])));
    if let Some(vsearch) = &vsearch {
        group.bench_function("vsearch", |b| b.iter(|| run(vsearch, &vargs)));
    }
    group.finish();
}

criterion_group!(benches, bench_binary);
criterion_main!(benches);
