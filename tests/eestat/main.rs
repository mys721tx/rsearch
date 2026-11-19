//! End-to-end tests of the eestat binary.
//!
//! Expected tables are those of `vsearch --fastq_eestats --output`,
//! tab-separated, a line for each position. A line gives the reads that
//! reach the position. It then gives the minimum, lower quartile, median,
//! mean, upper quartile and maximum of the quality scores, of their error
//! probabilities and of the expected errors up to that position. Scores
//! below zero count as zero. `--ascii`, `--qmin` and `--qmax` are vsearch's
//! `--fastq_ascii`, `--fastq_qmin` and `--fastq_qmax`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs eestat with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_eestat"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("eestat should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs eestat and returns its table, which must succeed.
fn eestat(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "eestat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
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
    std::env::temp_dir().join(format!("rsearch-eestat-{}-{name}", std::process::id()))
}

const HEADER: &str = "Pos Recs PctRecs \
                      Min_Q Low_Q Med_Q Mean_Q Hi_Q Max_Q \
                      Min_Pe Low_Pe Med_Pe Mean_Pe Hi_Pe Max_Pe \
                      Min_EE Low_EE Med_EE Mean_EE Hi_EE Max_EE";

/// The table of `rows`, written with single spaces between columns, under
/// the header; columns are separated by tabs.
fn table(rows: &str) -> String {
    std::iter::once(HEADER)
        .chain(rows.lines().map(str::trim).filter(|row| !row.is_empty()))
        .map(|row| row.split(' ').collect::<Vec<_>>().join("\t") + "\n")
        .collect()
}

#[test]
fn reports_each_position() {
    let input = fastq(&[
        ("a", "ACGTACGTAC", "IIIII5555#"),
        ("b", "ACGTAC", "III+++"),
        ("c", "AC", "!I"),
        ("d", "ACGTACGTAC", "IIIIIIIIII"),
    ]);
    let expected = table(
        "
        1 4 100.0 0.0 0.0 40.0 30.0 40.0 40.0 0.0001 0.0001 0.0001 0.25 0.0001 1 0.00 0.00 0.00 0.25 0.00 1.00
        2 4 100.0 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.00 0.00 0.00 0.25 0.00 1.00
        3 3 75.0 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.00 0.00 0.00 0.00 0.00 0.00
        4 3 75.0 10.0 10.0 40.0 30.0 40.0 40.0 0.0001 0.0001 0.0001 0.033 0.1 0.1 0.00 0.00 0.00 0.03 0.10 0.10
        5 3 75.0 10.0 10.0 40.0 30.0 40.0 40.0 0.0001 0.0001 0.0001 0.033 0.1 0.1 0.00 0.00 0.00 0.07 0.20 0.20
        6 3 75.0 10.0 10.0 20.0 23.3 40.0 40.0 0.0001 0.0001 0.01 0.037 0.1 0.1 0.00 0.00 0.01 0.10 0.30 0.30
        7 2 50.0 20.0 20.0 20.0 30.0 40.0 40.0 0.0001 0.0001 0.0001 0.005 0.01 0.01 0.00 0.00 0.00 0.01 0.02 0.02
        8 2 50.0 20.0 20.0 20.0 30.0 40.0 40.0 0.0001 0.0001 0.0001 0.005 0.01 0.01 0.00 0.00 0.00 0.02 0.03 0.03
        9 2 50.0 20.0 20.0 20.0 30.0 40.0 40.0 0.0001 0.0001 0.0001 0.005 0.01 0.01 0.00 0.00 0.00 0.02 0.04 0.04
        10 2 50.0 2.0 2.0 2.0 21.0 40.0 40.0 0.0001 0.0001 0.0001 0.32 0.63 0.63 0.00 0.00 0.00 0.34 0.67 0.67
        ",
    );
    assert_eq!(eestat(&[], &input), expected);
}

#[test]
fn bins_expected_errors_by_thousandths() {
    // quartiles of expected errors are the midpoints of bins of 0.001,
    // the mean is exact
    let input = fastq(&[
        ("a", "ACGTACGTACGTAC", "++++++++++++++"),
        ("b", "ACGTACGTACGT", "5555555555++"),
        ("c", "ACGTACGTACGTACGTACGTA", "+++++IIIII+++++IIIII#"),
    ]);
    let expected = table(
        "
        1 3 100.0 10.0 10.0 10.0 13.3 20.0 20.0 0.01 0.01 0.1 0.07 0.1 0.1 0.01 0.01 0.10 0.07 0.10 0.10
        2 3 100.0 10.0 10.0 10.0 13.3 20.0 20.0 0.01 0.01 0.1 0.07 0.1 0.1 0.02 0.02 0.20 0.14 0.20 0.20
        3 3 100.0 10.0 10.0 10.0 13.3 20.0 20.0 0.01 0.01 0.1 0.07 0.1 0.1 0.03 0.03 0.30 0.21 0.30 0.30
        4 3 100.0 10.0 10.0 10.0 13.3 20.0 20.0 0.01 0.01 0.1 0.07 0.1 0.1 0.04 0.04 0.40 0.28 0.40 0.40
        5 3 100.0 10.0 10.0 10.0 13.3 20.0 20.0 0.01 0.01 0.1 0.07 0.1 0.1 0.05 0.05 0.50 0.35 0.50 0.50
        6 3 100.0 10.0 10.0 20.0 23.3 40.0 40.0 0.0001 0.0001 0.01 0.037 0.1 0.1 0.06 0.06 0.50 0.39 0.60 0.60
        7 3 100.0 10.0 10.0 20.0 23.3 40.0 40.0 0.0001 0.0001 0.01 0.037 0.1 0.1 0.07 0.07 0.50 0.42 0.70 0.70
        8 3 100.0 10.0 10.0 20.0 23.3 40.0 40.0 0.0001 0.0001 0.01 0.037 0.1 0.1 0.08 0.08 0.50 0.46 0.80 0.80
        9 3 100.0 10.0 10.0 20.0 23.3 40.0 40.0 0.0001 0.0001 0.01 0.037 0.1 0.1 0.09 0.09 0.50 0.50 0.90 0.90
        10 3 100.0 10.0 10.0 20.0 23.3 40.0 40.0 0.0001 0.0001 0.01 0.037 0.1 0.1 0.10 0.10 0.50 0.53 1.00 1.00
        11 3 100.0 10.0 10.0 10.0 10.0 10.0 10.0 0.1 0.1 0.1 0.1 0.1 0.1 0.20 0.20 0.60 0.63 1.10 1.10
        12 3 100.0 10.0 10.0 10.0 10.0 10.0 10.0 0.1 0.1 0.1 0.1 0.1 0.1 0.30 0.30 0.70 0.73 1.20 1.20
        13 2 66.7 10.0 10.0 10.0 10.0 10.0 10.0 0.1 0.1 0.1 0.1 0.1 0.1 0.80 0.80 0.80 1.05 1.30 1.30
        14 2 66.7 10.0 10.0 10.0 10.0 10.0 10.0 0.1 0.1 0.1 0.1 0.1 0.1 0.90 0.90 0.90 1.15 1.40 1.40
        15 1 33.3 10.0 10.0 10.0 10.0 10.0 10.0 0.1 0.1 0.1 0.1 0.1 0.1 1.00 1.00 1.00 1.00 1.00 1.00
        16 1 33.3 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 1.00 1.00 1.00 1.00 1.00 1.00
        17 1 33.3 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 1.00 1.00 1.00 1.00 1.00 1.00
        18 1 33.3 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 1.00 1.00 1.00 1.00 1.00 1.00
        19 1 33.3 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 1.00 1.00 1.00 1.00 1.00 1.00
        20 1 33.3 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 1.00 1.00 1.00 1.00 1.00 1.00
        21 1 33.3 2.0 2.0 2.0 2.0 2.0 2.0 0.63 0.63 0.63 0.63 0.63 0.63 1.63 1.63 1.63 1.63 1.63 1.63
        ",
    );
    assert_eq!(eestat(&[], &input), expected);
}

#[test]
fn writes_probabilities_with_two_significant_digits() {
    // C's %.2g: plain down to 0.0001, then with an exponent
    let input = fastq(&[("a", "ACGTAC", "?5I!sK")]);
    let expected = table(
        "
        1 1 100.0 30.0 30.0 30.0 30.0 30.0 30.0 0.001 0.001 0.001 0.001 0.001 0.001 0.00 0.00 0.00 0.00 0.00 0.00
        2 1 100.0 20.0 20.0 20.0 20.0 20.0 20.0 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01
        3 1 100.0 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.01 0.01 0.01 0.01 0.01 0.01
        4 1 100.0 0.0 0.0 0.0 0.0 0.0 0.0 1 1 1 1 1 1 1.01 1.01 1.01 1.01 1.01 1.01
        5 1 100.0 82.0 82.0 82.0 82.0 82.0 82.0 6.3e-09 6.3e-09 6.3e-09 6.3e-09 6.3e-09 6.3e-09 1.01 1.01 1.01 1.01 1.01 1.01
        6 1 100.0 42.0 42.0 42.0 42.0 42.0 42.0 6.3e-05 6.3e-05 6.3e-05 6.3e-05 6.3e-05 6.3e-05 1.01 1.01 1.01 1.01 1.01 1.01
        ",
    );
    assert_eq!(eestat(&[], &input), expected);
}

#[test]
fn reports_empty_input_and_reads() {
    assert_eq!(eestat(&[], ""), table(""));
    assert_eq!(eestat(&[], &fastq(&[("e", "", "")])), table(""));
    // an empty read counts in the share of reads reaching each position
    let input = fastq(&[("e", "", ""), ("a", "A", "I")]);
    let expected = table(
        "1 1 50.0 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.00 0.00 0.00 0.00 0.00 0.00",
    );
    assert_eq!(eestat(&[], &input), expected);
}

#[test]
fn counts_negative_scores_as_zero() {
    let input = fastq(&[("a", "ACGT", ";@hh")]);
    let expected = table(
        "
        1 1 100.0 0.0 0.0 0.0 0.0 0.0 0.0 1 1 1 1 1 1 1.00 1.00 1.00 1.00 1.00 1.00
        2 1 100.0 0.0 0.0 0.0 0.0 0.0 0.0 1 1 1 1 1 1 2.00 2.00 2.00 2.00 2.00 2.00
        3 1 100.0 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 2.00 2.00 2.00 2.00 2.00 2.00
        4 1 100.0 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 2.00 2.00 2.00 2.00 2.00 2.00
        ",
    );
    assert_eq!(eestat(&["--ascii", "64", "--qmin", "-5"], &input), expected);
}

#[test]
fn checks_the_range_of_scores() {
    // scores up to 126 - offset are accepted by default
    let input = fastq(&[("a", "AC", "IK")]);
    assert!(eestat(&[], &input).contains("\t42.0\t"));
    // '!' is -31 at offset 64
    let low = fastq(&[("a", "AC", "I!")]);
    for (args, input, message) in [
        (&["--qmax", "41"][..], &input, "above"),
        (&["--qmin", "41"], &input, "below"),
        (&["--ascii", "64"], &low, "below"),
    ] {
        let output = run(args, input);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{args:?}: {stderr}");
    }
}

#[test]
fn ignores_abundances() {
    // vsearch --fastq_eestats reads no abundance and checks none
    let output = run(&[], &fastq(&[("a;size=0", "AC", "II")]));
    assert!(output.status.success());
}

#[test]
fn rejects_fasta() {
    let output = run(&[], ">a\nACGT\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("FASTQ file expected"), "{stderr}");
}

#[test]
fn rejects_invalid_input() {
    let output = run(&[], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let output = run(&[], "@a\nAC\n+\nI\x1f\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--unknown"][..],
        &["--ascii", "50"],
        &["--qmin", "10", "--qmax", "5"],
        &["--qmin", "-1"],
        &["--ascii", "64", "--qmax", "63"],
    ] {
        let output = run(args, &fastq(&[("a", "AC", "II")]));
        assert_eq!(output.status.code(), Some(2), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.tsv"));
    std::fs::write(&input, fastq(&[("a", "A", "I")])).unwrap();
    let stdout = eestat(
        &[
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
    let expected = table(
        "1 1 100.0 40.0 40.0 40.0 40.0 40.0 40.0 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.00 0.00 0.00 0.00 0.00 0.00",
    );
    assert_eq!(written, expected);
}
