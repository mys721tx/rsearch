//! End-to-end tests of the shufl binary.
//!
//! Expected outputs are those of `vsearch --shuffle --notrunclabels
//! --minseqlength 0 --fasta_width 0` with `--randseed` (`--seed`) on the
//! same input. A given seed gives vsearch's order. Without a seed, or with
//! 0, the seed comes from the system, and the order is checked only to be a
//! permutation of the input. vsearch writes only FASTA; FASTQ reads keep
//! their scores, and are checked by writing them as FASTA with conv.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs a binary with `args`, feeding `input` on stdin.
fn run_bin(bin: &str, args: &[&str], input: &str) -> Output {
    let mut child = Command::new(bin)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("the binary should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs shufl with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_shufl"), args, input)
}

/// Returns the standard output of a run, which must have succeeded.
fn stdout_of(name: &str, output: Output) -> String {
    assert!(
        output.status.success(),
        "{name} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

/// Runs shufl and returns its standard output, which must succeed.
fn shufl(args: &[&str], input: &str) -> String {
    stdout_of("shufl", run(args, input))
}

/// Builds FASTA text from (label, sequence) pairs.
fn fasta(records: &[(&str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq)| format!(">{id}\n{seq}\n"))
        .collect()
}

/// Builds FASTQ text from (label, sequence, quality) triples.
fn fastq(records: &[(&str, &str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq, qual)| format!("@{id}\n{seq}\n+\n{qual}\n"))
        .collect()
}

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-shufl-{}-{name}", std::process::id()))
}

/// Eight records, one with a description, one with an abundance and one
/// empty.
const RECORDS: [(&str, &str); 8] = [
    ("r0", "A"),
    ("r1 d", "CC"),
    ("r2;size=3", "GGG"),
    ("r3", "T"),
    ("r4", ""),
    ("r5", "AC"),
    ("r6", "AG"),
    ("r7", "AT"),
];

fn records() -> String {
    fasta(&RECORDS)
}

/// The records in the order of their indices.
fn ordered(indices: [usize; 8]) -> String {
    let records: Vec<_> = indices.iter().map(|&i| RECORDS[i]).collect();
    fasta(&records)
}

#[test]
fn shuffles_as_vsearch_with_a_seed() {
    for (seed, order) in [
        ("1", [3, 4, 5, 7, 6, 2, 0, 1]),
        ("2", [3, 2, 0, 1, 6, 4, 5, 7]),
        ("7", [2, 1, 3, 5, 4, 0, 7, 6]),
        ("42", [2, 1, 5, 3, 0, 7, 4, 6]),
        ("9223372036854775807", [0, 1, 6, 5, 2, 7, 3, 4]),
    ] {
        assert_eq!(
            shufl(&["--seed", seed], &records()),
            ordered(order),
            "{seed}"
        );
    }
}

#[test]
fn takes_negative_seeds_as_vsearch() {
    // -1 is the seed 2^64 - 1
    let expected = ordered([4, 6, 1, 3, 2, 7, 5, 0]);
    assert_eq!(shufl(&["--seed", "-1"], &records()), expected);
}

#[test]
fn shuffles_many_records_as_vsearch() {
    let records: Vec<(String, String)> = (0..30)
        .map(|i| (format!("s{i}"), "ACGT"[i % 4..i % 4 + 1].repeat(i % 5 + 1)))
        .collect();
    let input: String = records
        .iter()
        .map(|(l, s)| format!(">{l}\n{s}\n"))
        .collect();
    for (seed, order) in [
        (
            "1",
            [
                8, 23, 7, 20, 26, 19, 25, 17, 15, 21, 2, 5, 18, 29, 16, 6, 28, 14, 10, 24, 13, 27,
                1, 11, 22, 9, 0, 12, 3, 4,
            ],
        ),
        (
            "123456789",
            [
                15, 11, 20, 17, 27, 14, 2, 26, 19, 21, 25, 24, 9, 23, 5, 13, 6, 29, 16, 4, 18, 12,
                28, 1, 8, 22, 0, 3, 7, 10,
            ],
        ),
    ] {
        let expected: String = order
            .iter()
            .map(|&i| format!(">{}\n{}\n", records[i].0, records[i].1))
            .collect();
        assert_eq!(shufl(&["--seed", seed], &input), expected, "{seed}");
    }
}

/// Sorts the records of FASTA text, to compare permutations.
fn sorted(text: &str) -> Vec<String> {
    let mut records: Vec<String> = text
        .split('>')
        .filter(|r| !r.is_empty())
        .map(String::from)
        .collect();
    records.sort();
    records
}

#[test]
fn shuffles_with_a_system_seed_without_one() {
    for args in [&[][..], &["--seed", "0"]] {
        let output = shufl(args, &records());
        assert_eq!(sorted(&output), sorted(&records()), "{args:?}");
    }
}

#[test]
fn shuffles_fastq() {
    let input = fastq(&[
        ("q", "AC", "II"),
        ("r", "G", "5"),
        ("s", "TT", "#I"),
        ("t", "", ""),
    ]);
    let expected = fastq(&[
        ("q", "AC", "II"),
        ("s", "TT", "#I"),
        ("r", "G", "5"),
        ("t", "", ""),
    ]);
    let output = shufl(&["--seed", "42"], &input);
    assert_eq!(output, expected);
    // vsearch --shuffle writes FASTQ input as FASTA
    let conv = env!("CARGO_BIN_EXE_conv");
    let as_fasta = stdout_of(
        "conv",
        run_bin(conv, &["--fasta", "--fasta-width", "0"], &output),
    );
    let expected = fasta(&[("q", "AC"), ("s", "TT"), ("r", "G"), ("t", "")]);
    assert_eq!(as_fasta, expected);
}

#[test]
fn keeps_a_single_record() {
    let input = fasta(&[("a", "A")]);
    assert_eq!(shufl(&["--seed", "5"], &input), input);
}

#[test]
fn rejects_invalid_seeds() {
    // vsearch takes a signed 64-bit integer
    for seed in ["9223372036854775808", "18446744073709551615", "x", "1.5"] {
        let output = run(&["--seed", seed], &records());
        assert!(!output.status.success(), "{seed} should fail");
    }
}

#[test]
fn rejects_invalid_abundances_before_writing() {
    // every record is read, and checked, before any is written
    let output = run(&["--seed", "5"], &fasta(&[("b", "C"), ("a;size=0", "A")]));
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
}

#[test]
fn rejects_invalid_input_before_writing() {
    let output = run(&["--seed", "5"], ">a\nACGT\n>b\nAC-T\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character '-'"),
        "{stderr}"
    );
}

#[test]
fn unwraps_fasta_and_strips_invalid_characters() {
    let output = run(&["--seed", "5"], ">a\nAC\nG1T\n");
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a", "ACGT")])
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn accepts_empty_input() {
    assert_eq!(shufl(&["--seed", "5"], ""), "");
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, records()).unwrap();
    let stdout = shufl(
        &[
            "--seed",
            "1",
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(written, ordered([3, 4, 5, 7, 6, 2, 0, 1]));
}
