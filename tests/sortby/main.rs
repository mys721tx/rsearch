//! End-to-end tests of the sortby binary.
//!
//! Expected outputs are those of `vsearch --sortbylength` (`--length`) and
//! `vsearch --sortbysize` (`--size`) with `--notrunclabels --minseqlength 0
//! --fasta_width 0` on the same input. They sort by decreasing length, then
//! abundance, then whole header, then input order; or by decreasing
//! abundance, then whole header, then input order. Headers are written
//! unchanged. vsearch writes only FASTA; FASTQ reads keep their scores, and
//! are checked by writing them as FASTA with conv.

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
    let mut stdin = child.stdin.take().unwrap();
    let fed = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// Runs sortby with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_sortby"), args, input)
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

/// Runs sortby and returns its standard output, which must succeed.
fn sortby(args: &[&str], input: &str) -> String {
    stdout_of("sortby", run(args, input))
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
    std::env::temp_dir().join(format!("rsearch-sortby-{}-{name}", std::process::id()))
}

const LONG: &str = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

/// Reads with ties at every level: two with the same header and abundance,
/// two headers that differ only in case, an abundance a description hides,
/// and an empty sequence.
fn reads() -> String {
    fasta(&[
        ("b;size=2", "ACG"),
        ("a", "ACG"),
        ("c;size=5 x", "A"),
        ("d;size=3", "ACGTACGT"),
        ("e", ""),
        ("a", "TT"),
        ("A", "GGG"),
        ("b;size=2", "CCC"),
        ("x;size=9", LONG),
    ])
}

#[test]
fn sorts_by_length() {
    // then by abundance, then header (upper case first), then input order
    let expected = fasta(&[
        ("x;size=9", LONG),
        ("d;size=3", "ACGTACGT"),
        ("b;size=2", "ACG"),
        ("b;size=2", "CCC"),
        ("A", "GGG"),
        ("a", "ACG"),
        ("a", "TT"),
        ("c;size=5 x", "A"),
        ("e", ""),
    ]);
    assert_eq!(sortby(&["--length"], &reads()), expected);
}

#[test]
fn sorts_by_abundance() {
    // then header, then input order; a size followed by a description is
    // not read, and counts as 1
    let expected = fasta(&[
        ("x;size=9", LONG),
        ("d;size=3", "ACGTACGT"),
        ("b;size=2", "ACG"),
        ("b;size=2", "CCC"),
        ("A", "GGG"),
        ("a", "ACG"),
        ("a", "TT"),
        ("c;size=5 x", "A"),
        ("e", ""),
    ]);
    assert_eq!(sortby(&["--size"], &reads()), expected);
}

#[test]
fn sorts_abundances_beyond_32_bits() {
    let input = fasta(&[
        ("c;size=4294967295", "A"),
        ("b;size=4294967296", "ACG"),
        ("a;size=9223372036854775807", "AC"),
    ]);
    let expected = fasta(&[
        ("a;size=9223372036854775807", "AC"),
        ("b;size=4294967296", "ACG"),
        ("c;size=4294967295", "A"),
    ]);
    assert_eq!(sortby(&["--size"], &input), expected);
}

#[test]
fn sorts_fastq() {
    let input = fastq(&[
        ("q;size=2", "AC", "I#"),
        ("r", "ACGT", "II5I"),
        ("s;size=4", "AC", "55"),
        ("e", "", ""),
    ]);
    let expected = fastq(&[
        ("r", "ACGT", "II5I"),
        ("s;size=4", "AC", "55"),
        ("q;size=2", "AC", "I#"),
        ("e", "", ""),
    ]);
    assert_eq!(sortby(&["--length"], &input), expected);
    let expected = fastq(&[
        ("s;size=4", "AC", "55"),
        ("q;size=2", "AC", "I#"),
        ("e", "", ""),
        ("r", "ACGT", "II5I"),
    ]);
    assert_eq!(sortby(&["--size"], &input), expected);
}

#[test]
fn matches_vsearch_on_fastq_written_as_fasta() {
    // vsearch --sortbylength and --sortbysize read FASTQ and write FASTA
    let input = fastq(&[
        ("q;size=2", "AC", "I#"),
        ("r", "ACGT", "II5I"),
        ("s;size=4", "AC", "55"),
        ("e", "", ""),
    ]);
    let conv = env!("CARGO_BIN_EXE_conv");
    let as_fasta = |text: String| {
        stdout_of(
            "conv",
            run_bin(conv, &["--fasta", "--fasta-width", "0"], &text),
        )
    };
    let expected = fasta(&[
        ("r", "ACGT"),
        ("s;size=4", "AC"),
        ("q;size=2", "AC"),
        ("e", ""),
    ]);
    assert_eq!(as_fasta(sortby(&["--length"], &input)), expected);
    let expected = fasta(&[
        ("s;size=4", "AC"),
        ("q;size=2", "AC"),
        ("e", ""),
        ("r", "ACGT"),
    ]);
    assert_eq!(as_fasta(sortby(&["--size"], &input)), expected);
}

#[test]
fn unwraps_fasta_and_strips_invalid_characters() {
    let output = run(&["--length"], ">b\nAC1\n>a\nAC\nG\n");
    assert!(output.status.success());
    let expected = fasta(&[("a", "ACG"), ("b", "AC")]);
    assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn rejects_invalid_abundances_before_writing() {
    // every read is read, and checked, before any is written
    for (input, message) in [
        (fasta(&[("a;size=1", "AC"), ("z;size=0", "AC")]), "zero"),
        (
            fasta(&[("a;size=1", "AC"), ("z;size=9223372036854775808", "AC")]),
            "range",
        ),
    ] {
        for key in ["--length", "--size"] {
            let output = run(&[key], &input);
            assert_eq!(output.status.code(), Some(1), "{key}");
            assert!(output.stdout.is_empty(), "{key}");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(stderr.contains(message), "{key}: {stderr}");
        }
    }
}

#[test]
fn rejects_invalid_input_before_writing() {
    let output = run(&["--length"], ">a\nACGT\n>b\nAC-T\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character '-' on line 4 of FASTA file"),
        "{stderr}"
    );
}

#[test]
fn needs_exactly_one_key() {
    for args in [
        &[][..],
        &["--length", "--size"],
        &["--label"],
        &["--topn", "2"],
    ] {
        assert!(
            !run(args, &reads()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn accepts_empty_input() {
    assert_eq!(sortby(&["--length"], ""), "");
    assert_eq!(sortby(&["--size"], ""), "");
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, fasta(&[("a", "A"), ("b;size=2", "AC")])).unwrap();
    let stdout = sortby(
        &[
            "--size",
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
    assert_eq!(written, fasta(&[("b;size=2", "AC"), ("a", "A")]));
}

/// Runs a chain of rsearch binaries on `input`, each with its arguments,
/// and returns the last one's output.
fn chain(steps: &[(&str, &[&str])], input: &str) -> String {
    steps.iter().fold(input.to_string(), |text, (bin, args)| {
        stdout_of(bin, run_bin(bin, args, &text))
    })
}

#[test]
fn pipeline_sortbysize() {
    // vsearch --sortbysize in --notrunclabels --minsize 2 --topn 3 --output
    // out, whose default --minseqlength 1 drops the empty read
    let steps: [(&str, &[&str]); 4] = [
        (env!("CARGO_BIN_EXE_sqfilt"), &["--min-size", "2"]),
        (env!("CARGO_BIN_EXE_sortby"), &["--size"]),
        (env!("CARGO_BIN_EXE_topn"), &["3"]),
        (env!("CARGO_BIN_EXE_conv"), &["--fasta"]),
    ];
    let expected = fasta(&[
        ("x;size=9", LONG),
        ("d;size=3", "ACGTACGT"),
        ("b;size=2", "ACG"),
    ]);
    assert_eq!(chain(&steps, &reads()), expected);
}
