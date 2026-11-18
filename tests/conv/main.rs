//! End-to-end tests of the conv binary.
//!
//! Expected outputs are those of `vsearch --fastq_convert` (FASTQ input) and
//! `vsearch --fasta2fastq` (FASTA input) on the same input, with the
//! equivalent options.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs conv with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_conv"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("conv should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs conv and returns its standard output, which must succeed.
fn conv(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "conv failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
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
    std::env::temp_dir().join(format!("rsearch-conv-{}-{name}", std::process::id()))
}

#[test]
fn copies_phred33_unchanged() {
    // the whole header is kept, description included
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    assert_eq!(conv(&[], &input), input);
}

#[test]
fn writes_phred64() {
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    let expected = fastq(&[("a;size=3 desc", "ACGT", "@JT~")]);
    let args = ["--ascii-out", "64", "--qmax-out", "62"];
    assert_eq!(conv(&args, &input), expected);
}

#[test]
fn reads_phred64() {
    let input = fastq(&[("a", "ACGT", "@Jh~")]);
    let expected = fastq(&[("a", "ACGT", "!+I_")]);
    assert_eq!(conv(&["--ascii", "64"], &input), expected);
}

#[test]
fn limits_written_scores() {
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    let expected = fastq(&[("a;size=3 desc", "ACGT", "&+5I")]);
    let args = ["--qmin-out", "5", "--qmax-out", "40"];
    assert_eq!(conv(&args, &input), expected);
}

#[test]
fn converts_solexa_to_phred() {
    // Solexa -5, -4, -1, 0, 1, 2, 10, 40
    let input = fastq(&[("s", "ACGTACGT", ";<?@ACJh")]);
    let expected = fastq(&[("s", "ACGTACGT", "\"\"$$%&+I")]);
    assert_eq!(conv(&["--ascii", "64", "--solexa"], &input), expected);
}

#[test]
fn relabels() {
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    let expected = fastq(&[("R1;size=1", "ACGT", "!+5~")]);
    assert_eq!(conv(&["--relabel", "R", "--sizeout"], &input), expected);
}

#[test]
fn writes_fasta_as_fastq() {
    let input = fasta(&[("a;size=3 desc", "ACGT"), ("b", "AC")]);
    let expected = fastq(&[("a;size=3 desc", "ACGT", "JJJJ"), ("b", "AC", "JJ")]);
    assert_eq!(conv(&[], &input), expected);
    // the highest output score, written with the output offset
    let expected = fastq(&[
        ("a;size=3 desc;size=1", "ACGT", "rrrr"),
        ("b;size=1", "AC", "rr"),
    ]);
    let args = ["--ascii-out", "64", "--qmax-out", "50", "--sizeout"];
    assert_eq!(conv(&args, &input), expected);
}

#[test]
fn stops_at_scores_out_of_range() {
    // the reads before are written
    let input = fastq(&[("a", "AC", "II"), ("b", "ACGT", "!+5~")]);
    let output = run(&["--qmax", "41"], &input);
    assert!(!output.status.success());
    assert_eq!(output.stdout, fastq(&[("a", "AC", "II")]).as_bytes());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("quality score 93 above the maximum 41 in read b"));
}

#[test]
fn rejects_invalid_encodings() {
    let input = fastq(&[("a", "ACGT", "IIII")]);
    for args in [
        &["--solexa"][..],
        &["--ascii", "40"],
        &["--qmin", "10", "--qmax", "5"],
        &["--ascii", "64", "--qmin", "-32"],
        &["--qmax", "94"],
        &["--qmin-out", "10", "--qmax-out", "5"],
        &["--qmin-out", "-1"],
        &["--ascii-out", "64", "--qmax-out", "63"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn rejects_input_score_options_with_fasta() {
    let input = fasta(&[("a", "ACGT")]);
    for args in [&["--ascii", "64"][..], &["--qmin", "0"], &["--solexa"]] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, fastq(&[("a", "ACGT", "@Jh~")])).unwrap();
    let stdout = conv(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--ascii",
            "64",
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(written, fastq(&[("a", "ACGT", "!+I_")]));
}

#[test]
fn accepts_empty_input() {
    assert_eq!(conv(&[], ""), "");
}
