//! End-to-end tests of the sqtrim binary.
//!
//! Expected outputs are those of `vsearch --fastx_filter` with the
//! equivalent trimming options. Its `--fastqout` (or `--fastaout
//! --fasta_width 0`) and its `--fastqout_discarded` (or
//! `--fastaout_discarded`) are taken together, in input order. sqtrim never
//! discards. The reads vsearch discards because trimming emptied them are
//! therefore written, empty, as vsearch writes them to its discarded
//! output. Lengths are truncated as with `--fastq_trunclen_keep`;
//! `--fastq_trunclen`, which also discards shorter reads, is `sqtrim --len`
//! then `sqfilt --min-len`.

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

/// Runs sqtrim with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_sqtrim"), args, input)
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

/// Runs sqtrim and returns its standard output, which must succeed.
fn sqtrim(args: &[&str], input: &str) -> String {
    stdout_of("sqtrim", run(args, input))
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
    std::env::temp_dir().join(format!("rsearch-sqtrim-{}-{name}", std::process::id()))
}

/// Four reads: a long one with low scores near its end, a short one of high
/// scores, one of two bases, and one that starts with two scores of 0.
const LABELS: [&str; 4] = ["a desc", "b", "c", "e"];

fn reads() -> String {
    fastq(&[
        ("a desc", "ACGTACGTAC", "IIIII5+#II"),
        ("b", "ACGTA", "IIIII"),
        ("c", "AC", "#I"),
        ("e", "ACGT", "!!II"),
    ])
}

/// The four reads trimmed to the given sequences and scores.
fn trimmed(reads: [(&str, &str); 4]) -> String {
    let records: Vec<_> = (0..4)
        .map(|i| (LABELS[i], reads[i].0, reads[i].1))
        .collect();
    fastq(&records)
}

#[test]
fn copies_reads_unchanged() {
    assert_eq!(sqtrim(&[], &reads()), reads());
}

#[test]
fn strips_the_left_end() {
    let expected = trimmed([
        ("GTACGTAC", "III5+#II"),
        ("GTA", "III"),
        ("", ""),
        ("GT", "II"),
    ]);
    assert_eq!(sqtrim(&["--strip-left", "2"], &reads()), expected);
}

#[test]
fn strips_the_right_end() {
    let expected = trimmed([("ACGTACG", "IIIII5+"), ("AC", "II"), ("", ""), ("A", "!")]);
    assert_eq!(sqtrim(&["--strip-right", "3"], &reads()), expected);
}

#[test]
fn strips_both_ends() {
    // reads shorter than both strips together are emptied
    let expected = trimmed([("TACG", "II5+"), ("", ""), ("", ""), ("", "")]);
    let args = ["--strip-left", "3", "--strip-right", "3"];
    assert_eq!(sqtrim(&args, &reads()), expected);
}

#[test]
fn truncates_to_a_length() {
    // vsearch --fastq_trunclen_keep: shorter reads are kept as they are
    let expected = trimmed([
        ("ACGT", "IIII"),
        ("ACGT", "IIII"),
        ("AC", "#I"),
        ("ACGT", "!!II"),
    ]);
    assert_eq!(sqtrim(&["--len", "4"], &reads()), expected);
}

#[test]
fn truncates_at_a_low_score() {
    // before the first base with a score at or below the limit
    let expected = trimmed([("ACGTAC", "IIIII5"), ("ACGTA", "IIIII"), ("", ""), ("", "")]);
    assert_eq!(sqtrim(&["--qual", "10"], &reads()), expected);
    let expected = trimmed([
        ("ACGTACG", "IIIII5+"),
        ("ACGTA", "IIIII"),
        ("", ""),
        ("", ""),
    ]);
    assert_eq!(sqtrim(&["--qual", "2"], &reads()), expected);
}

#[test]
fn truncates_at_an_expected_error_limit() {
    // before the base that takes the expected errors above the limit
    let expected = trimmed([
        ("ACGTACG", "IIIII5+"),
        ("ACGTA", "IIIII"),
        ("", ""),
        ("", ""),
    ]);
    assert_eq!(sqtrim(&["--max-ee", "0.5"], &reads()), expected);
    // every base has some probability of error
    let expected = trimmed([("", ""), ("", ""), ("", ""), ("", "")]);
    assert_eq!(sqtrim(&["--max-ee", "0"], &reads()), expected);
}

#[test]
fn truncates_at_an_expected_error_rate() {
    let expected = trimmed([
        ("ACGTACG", "IIIII5+"),
        ("ACGTA", "IIIII"),
        ("", ""),
        ("", ""),
    ]);
    assert_eq!(sqtrim(&["--max-ee-rate", "0.05"], &reads()), expected);
}

#[test]
fn strips_then_truncates() {
    // ends first, then the length, then base by base
    let expected = trimmed([("CGTAC", "IIII5"), ("CGTA", "IIII"), ("C", "I"), ("", "")]);
    let args = ["--strip-left", "1", "--len", "6", "--qual", "10"];
    assert_eq!(sqtrim(&args, &reads()), expected);
}

#[test]
fn checks_scores_of_the_bases_it_reads() {
    // as vsearch: scores out of range fail, but only those of the bases
    // left after stripping and truncating to a length, and before a
    // truncating base
    for (args, message) in [
        (
            &["--qual", "10", "--qmax", "39"][..],
            "quality score 40 above the maximum 39",
        ),
        (&["--ascii", "64"], "quality score -11 below the minimum 0"),
    ] {
        let output = run(args, &reads());
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{args:?}: {stderr}");
    }
    let input = fastq(&[("q", "ACGT", "~III")]);
    assert_eq!(run(&["--qmax", "41"], &input).status.code(), Some(1));
    let args = ["--qmax", "41", "--strip-left", "1"];
    assert_eq!(sqtrim(&args, &input), fastq(&[("q", "CGT", "III")]));
    let input = fastq(&[("q", "ACGT", "II#~")]);
    for (args, expected) in [
        (&["--qmax", "41", "--qual", "2"][..], ("AC", "II")),
        (&["--qmax", "41", "--len", "3"], ("ACG", "II#")),
        (&["--qmax", "41", "--max-ee", "0.1"], ("AC", "II")),
    ] {
        let expected = fastq(&[("q", expected.0, expected.1)]);
        assert_eq!(sqtrim(args, &input), expected, "{args:?}");
    }
}

#[test]
fn trims_fasta_by_position() {
    let input = fasta(&[("a x", "ACGTACGTAC"), ("b", "ACG")]);
    let args = ["--strip-left", "2", "--strip-right", "1"];
    let expected = fasta(&[("a x", "GTACGTA"), ("b", "")]);
    assert_eq!(sqtrim(&args, &input), expected);
    let expected = fasta(&[("a x", "ACGT"), ("b", "ACG")]);
    assert_eq!(sqtrim(&["--len", "4"], &input), expected);
}

#[test]
fn refuses_score_options_with_fasta() {
    let input = fasta(&[("a", "ACGT")]);
    for args in [
        &["--qual", "2"][..],
        &["--max-ee", "1"],
        &["--max-ee-rate", "0.1"],
        &["--ascii", "64"],
        &["--qmin", "0"],
        &["--qmax", "41"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn rejects_invalid_arguments() {
    let input = reads();
    for args in [
        &["--len", "0"][..],
        &["--len", "x"],
        &["--strip-left", "-1"],
        &["--max-ee", "-1"],
        &["--max-ee-rate", "-0.1"],
        &["--qual", "-1"],
        &["--ascii", "40"],
        &["--qmin", "10", "--qmax", "5"],
        &["--ascii", "64", "--qmax", "63"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn keeps_empty_reads() {
    let input = fastq(&[("e", "", ""), ("a", "ACGT", "IIII")]);
    let expected = fastq(&[("e", "", ""), ("a", "AC", "II")]);
    assert_eq!(sqtrim(&["--strip-right", "2"], &input), expected);
}

#[test]
fn accepts_empty_input() {
    assert_eq!(sqtrim(&[], ""), "");
    assert_eq!(sqtrim(&["--qual", "2"], ""), "");
}

#[test]
fn rejects_invalid_input() {
    let output = run(&["--len", "2"], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character 'I'"),
        "{stderr}"
    );
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, reads()).unwrap();
    let stdout = sqtrim(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--len",
            "4",
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(
        written,
        trimmed([
            ("ACGT", "IIII"),
            ("ACGT", "IIII"),
            ("AC", "#I"),
            ("ACGT", "!!II"),
        ])
    );
}

/// Runs sqtrim with `first` on `input`, then relab with `second` on its
/// output, and returns relab's output.
fn then_relab(first: &[&str], second: &[&str], input: &str) -> String {
    let trimmed = sqtrim(first, input);
    stdout_of(
        "relab",
        run_bin(env!("CARGO_BIN_EXE_relab"), second, &trimmed),
    )
}

#[test]
fn pipeline_expected_errors_of_trimmed_reads() {
    // vsearch --fastx_filter in --fastq_truncee 0.5 --eeout: the expected
    // errors of the read as written
    let input = fastq(&[
        ("a desc", "ACGTACGTAC", "IIIII5+#II"),
        ("b", "ACGTA", "IIIII"),
    ]);
    let expected = fastq(&[
        ("a desc;ee=0.1105", "ACGTACG", "IIIII5+"),
        ("b;ee=0.0005000", "ACGTA", "IIIII"),
    ]);
    assert_eq!(
        then_relab(&["--max-ee", "0.5"], &["--ee"], &input),
        expected
    );
    // vsearch --fastx_filter in --fastq_stripleft 1 --fastq_truncqual 10
    // --eeout --relabel T
    let expected = fastq(&[
        ("T1;ee=0.01040", "CGTAC", "IIII5"),
        ("T2;ee=0.0004000", "CGTA", "IIII"),
    ]);
    let args = ["--strip-left", "1", "--qual", "10"];
    assert_eq!(
        then_relab(&args, &["--ee", "--relabel", "T"], &input),
        expected
    );
}

#[test]
fn rejects_invalid_abundances() {
    // vsearch --fastx_filter reads the abundance of every read, and writes
    // the reads before the one that fails, trimmed: here kept and
    // discarded together
    for (bad, message) in [
        ("z;size=0", "zero"),
        ("z;size=9223372036854775808", "range"),
    ] {
        let input = fastq(&[
            ("a;size=1", "ACGT", "IIII"),
            ("b", "AC", "II"),
            (bad, "ACGT", "IIII"),
            ("c", "ACGT", "IIII"),
        ]);
        let output = run(&["--strip-left", "2"], &input);
        assert_eq!(output.status.code(), Some(1), "{bad}");
        let expected = fastq(&[("a;size=1", "GT", "II"), ("b", "", "")]);
        assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{bad}: {stderr}");
    }
    let input = fasta(&[("a;size=1", "ACGT"), ("z;size=0", "AC")]);
    let output = run(&["--len", "2"], &input);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a;size=1", "AC")])
    );
}
