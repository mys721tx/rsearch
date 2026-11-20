//! Times two builds of the same tools on the benchmark inputs, alternating
//! runs so that load on the machine affects both alike (`rscompare time`).
//!
//! The inputs are the files the criterion benchmarks write to `target/tmp`.
//! Run each benchmark once first (`cargo bench -- --test`).

use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Instant;

/// A timed command: its name, the tool, its input (under `target/tmp`) and
/// its other arguments.
const CASES: &[(&str, &str, &str, &[&str])] = &[
    ("derep fasta", "derep", "derep-bench.fa", &[]),
    ("derep fastq", "derep", "derep-bench.fq", &[]),
    (
        "sqtrim all",
        "sqtrim",
        "sqtrim-bench.fq",
        &["--strip-left", "5", "--len", "200", "--qual", "5"],
    ),
    (
        "sqfilt maxee",
        "sqfilt",
        "sqfilt-bench.fq",
        &["--max-ee", "1"],
    ),
    (
        "sqfilt fasta",
        "sqfilt",
        "sqfilt-bench.fa",
        &["--min-len", "100"],
    ),
    (
        "relab prefix",
        "relab",
        "relab-bench.fq",
        &["--relabel", "R", "--sizeout"],
    ),
    (
        "conv phred64",
        "conv",
        "conv-bench.fq",
        &["--ascii-out", "64", "--qmax-out", "62"],
    ),
    ("conv fasta", "conv", "conv-bench.fa", &[]),
    ("conv sff", "conv", "conv-bench.sff", &[]),
];

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn once(exe: &Path, input: &Path, args: &[&str]) -> std::io::Result<f64> {
    let start = Instant::now();
    let status = Command::new(exe)
        .arg("-i")
        .arg(input)
        .args(["-o", "/dev/null"])
        .args(args)
        .stderr(Stdio::null())
        .status()?;
    if !status.success() {
        return Err(std::io::Error::other(format!("{} failed", exe.display())));
    }
    Ok(start.elapsed().as_secs_f64())
}

/// Times each case `runs` times with each build (`baseline` and `new` hold
/// the tools), alternating which goes first, and prints the medians.
pub fn run(baseline: &Path, new: &Path, inputs: &Path, runs: usize) -> std::io::Result<()> {
    for (name, tool, input, args) in CASES {
        let input = inputs.join(input);
        let (mut base, mut fresh) = (Vec::new(), Vec::new());
        for i in 0..runs {
            if i % 2 == 1 {
                base.push(once(&baseline.join(tool), &input, args)?);
                fresh.push(once(&new.join(tool), &input, args)?);
            } else {
                fresh.push(once(&new.join(tool), &input, args)?);
                base.push(once(&baseline.join(tool), &input, args)?);
            }
        }
        let (b, n) = (median(base) * 1000.0, median(fresh) * 1000.0);
        println!(
            "{name:13} base {b:6.1} ms  new {n:6.1} ms  {:+5.1}%",
            100.0 * (n - b) / b
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn median() {
        assert_eq!(super::median(vec![3.0, 1.0, 2.0]), 2.0);
        assert_eq!(super::median(vec![4.0, 1.0, 2.0, 3.0]), 2.5);
    }
}
