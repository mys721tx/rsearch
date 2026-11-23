//! End-to-end tests of the eegrid binary.
//!
//! Expected grids are those of `vsearch --fastq_eestats2 --output`. A grid
//! gives the number of reads and their longest and average length. Then,
//! for each length cutoff, it gives the reads at least that long whose
//! expected errors up to it are within each expected-error cutoff. The
//! length cutoffs run from `--min-len` (default 50) by `--len-step`
//! (default 50) up to the longest read or `--max-len`, vsearch's
//! `--length_cutoffs`. `--max-ee` lists the expected-error cutoffs (default
//! 0.5,1,2), vsearch's `--ee_cutoffs`. `--ascii`, `--qmin` and `--qmax` are
//! vsearch's `--fastq_ascii`, `--fastq_qmin` and `--fastq_qmax`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs eegrid with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_eegrid"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("eegrid should start");
    let mut stdin = child.stdin.take().unwrap();
    let fed = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// Runs eegrid and returns its grid, which must succeed.
fn eegrid(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "eegrid failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-eegrid-{}-{name}", std::process::id()))
}

/// Reads of 10, 6, 2 and 10 bases.
fn reads() -> String {
    fastq(&[
        ("a", "ACGTACGTAC", "IIIII5555#"),
        ("b", "ACGTAC", "III+++"),
        ("c", "AC", "!I"),
        ("d", "ACGTACGTAC", "IIIIIIIIII"),
    ])
}

/// The header of a grid with the default expected-error cutoffs.
const DEFAULT_HEADER: &str = "\
Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
";

#[test]
fn starts_at_fifty_by_default() {
    // one row even when every read is shorter than the first cutoff
    let expected = format!(
        "4 reads, max len 10, avg 7.0\n\n{DEFAULT_HEADER}    \
         50          0(  0.0%)          0(  0.0%)          0(  0.0%)\n"
    );
    assert_eq!(eegrid(&[], &reads()), expected);
    let expected = format!(
        "4 reads, max len 10, avg 7.0\n\n{DEFAULT_HEADER}    \
         12          0(  0.0%)          0(  0.0%)          0(  0.0%)\n"
    );
    assert_eq!(
        eegrid(&["--min-len", "12", "--len-step", "1"], &reads()),
        expected
    );
}

#[test]
fn counts_reads_within_each_cutoff() {
    // the length cutoffs stop at the longest read
    let expected = "\
4 reads, max len 10, avg 7.0

Length         MaxEE 0.01         MaxEE 0.25         MaxEE 1.00
------   ----------------   ----------------   ----------------
     2          3( 75.0%)          3( 75.0%)          3( 75.0%)
     5          2( 50.0%)          3( 75.0%)          3( 75.0%)
     8          1( 25.0%)          2( 50.0%)          2( 50.0%)
";
    let args = [
        "--min-len",
        "2",
        "--len-step",
        "3",
        "--max-ee",
        "0.01,0.25,1",
    ];
    assert_eq!(eegrid(&args, &reads()), expected);
}

#[test]
fn stops_at_the_last_length_cutoff() {
    let expected = format!(
        "4 reads, max len 10, avg 7.0\n\n{DEFAULT_HEADER}     \
         2          3( 75.0%)          3( 75.0%)          4(100.0%)\n     \
         4          3( 75.0%)          3( 75.0%)          3( 75.0%)\n     \
         6          3( 75.0%)          3( 75.0%)          3( 75.0%)\n"
    );
    let args = ["--min-len", "2", "--max-len", "7", "--len-step", "2"];
    assert_eq!(eegrid(&args, &reads()), expected);
}

#[test]
fn keeps_the_order_of_expected_error_cutoffs() {
    let input = fastq(&[
        ("a", "ACGTACGTACGTAC", "++++++++++++++"),
        ("b", "ACGTACGTACGT", "5555555555++"),
        ("c", "ACGTACGTACGTACGTACGTA", "+++++IIIII+++++IIIII#"),
    ]);
    let expected = "\
3 reads, max len 21, avg 15.7

Length         MaxEE 1.00         MaxEE 0.50
------   ----------------   ----------------
     1          3(100.0%)          3(100.0%)
     5          3(100.0%)          3(100.0%)
     9          3(100.0%)          1( 33.3%)
    13          1( 33.3%)          0(  0.0%)
";
    let args = [
        "--min-len",
        "1",
        "--max-len",
        "15",
        "--len-step",
        "4",
        "--max-ee",
        "1,0.5",
    ];
    assert_eq!(eegrid(&args, &input), expected);
}

#[test]
fn reports_empty_input_and_reads() {
    assert_eq!(eegrid(&[], ""), format!("0 reads\n\n{DEFAULT_HEADER}"));
    // the longest read decides the rows, and empty reads give none
    assert_eq!(
        eegrid(&[], &fastq(&[("e", "", "")])),
        format!("1 reads, max len 0, avg 0.0\n\n{DEFAULT_HEADER}")
    );
    let input = fastq(&[("e", "", ""), ("a", "A", "I")]);
    let expected = format!(
        "2 reads, max len 1, avg 0.5\n\n{DEFAULT_HEADER}     \
         1          1( 50.0%)          1( 50.0%)          1( 50.0%)\n"
    );
    assert_eq!(
        eegrid(&["--min-len", "1", "--len-step", "1"], &input),
        expected
    );
}

#[test]
fn counts_negative_scores_as_certain_errors() {
    let input = fastq(&[("a", "ACGT", ";@hh")]);
    let expected = "\
1 reads, max len 4, avg 4.0

Length         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------
     1          1(100.0%)          1(100.0%)
     2          0(  0.0%)          1(100.0%)
     3          0(  0.0%)          0(  0.0%)
     4          0(  0.0%)          0(  0.0%)
";
    let args = [
        "--min-len",
        "1",
        "--len-step",
        "1",
        "--max-ee",
        "1,2",
        "--ascii",
        "64",
        "--qmin",
        "-5",
    ];
    assert_eq!(eegrid(&args, &input), expected);
}

#[test]
fn checks_the_range_of_scores() {
    // scores up to 126 - offset are accepted by default
    let input = fastq(&[("a", "AC", "IK")]);
    assert!(eegrid(&[], &input).starts_with("1 reads, max len 2"));
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
    // vsearch --fastq_eestats2 reads no abundance and checks none
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
    // vsearch: lengths of at least 1, the first no more than the last, a
    // step of at least 1, and expected errors above 0
    for args in [
        &["--unknown"][..],
        &["--min-len", "0"],
        &["--len-step", "0"],
        &["--min-len", "10", "--max-len", "9"],
        &["--max-ee", "0"],
        &["--max-ee", "1,-1"],
        &["--max-ee", "1,x"],
        &["--max-ee", "inf"],
        &["--ascii", "50"],
        &["--qmin", "10", "--qmax", "5"],
        &["--ascii", "64", "--qmax", "63"],
    ] {
        let output = run(args, &fastq(&[("a", "AC", "II")]));
        assert_eq!(output.status.code(), Some(2), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.txt"));
    std::fs::write(&input, reads()).unwrap();
    let stdout = eegrid(
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
    assert!(
        written.starts_with("4 reads, max len 10, avg 7.0\n\n"),
        "{written}"
    );
}
