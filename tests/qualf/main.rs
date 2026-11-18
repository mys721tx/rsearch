//! End-to-end tests of the qualf binary.
//!
//! Unless noted otherwise, expected outputs are those of
//! `vsearch --fastq_filter` with the corresponding options on the same
//! input, so that qualf keeps vsearch's filtering behaviour.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Reads covering the filters: high quality, a low-quality tail, two Ns
/// with an abundance, very low quality with an abundance hidden by a
/// description, a stale `ee=` annotation, and an empty read.
const INPUT: &str = "\
@r1\nACGTACGTAC\n+\nIIIIIIIIII\n\
@r2\nACGTACGTAC\n+\nIIIII+++++\n\
@r3;size=3\nACGTNNACGT\n+\n5555555555\n\
@r4;size=3 desc\nacgtn\n+\n#####\n\
@r5;ee=0.5;\nACG\n+\nIII\n\
@e\n\n+\n\n";

const R1: (&str, &str, &str) = ("r1", "ACGTACGTAC", "IIIIIIIIII");
const R2: (&str, &str, &str) = ("r2", "ACGTACGTAC", "IIIII+++++");
const R3: (&str, &str, &str) = ("r3;size=3", "ACGTNNACGT", "5555555555");
const R4: (&str, &str, &str) = ("r4;size=3 desc", "acgtn", "#####");
const R5: (&str, &str, &str) = ("r5;ee=0.5;", "ACG", "III");
const E: (&str, &str, &str) = ("e", "", "");

/// Runs qualf with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_qualf"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("qualf should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-qualf-{}-{name}", std::process::id()))
}

/// Runs qualf on `input` and returns the kept and the discarded reads.
fn qualf(args: &[&str], input: &str) -> (String, String) {
    let discard = temp_path(&format!("discard-{}", args.join("_")));
    let mut all_args = args.to_vec();
    all_args.extend(["--discard", discard.to_str().unwrap()]);
    let output = run(&all_args, input);
    assert!(
        output.status.success(),
        "qualf failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let discarded = std::fs::read_to_string(&discard).unwrap();
    std::fs::remove_file(discard).unwrap();
    (String::from_utf8(output.stdout).unwrap(), discarded)
}

/// Builds FASTQ text from (header, sequence, quality) triples.
fn fastq(records: &[(&str, &str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq, qual)| format!("@{id}\n{seq}\n+\n{qual}\n"))
        .collect()
}

#[track_caller]
fn check(args: &[&str], kept: &[(&str, &str, &str)], discarded: &[(&str, &str, &str)]) {
    let (out, dis) = qualf(args, INPUT);
    assert_eq!(out, fastq(kept), "kept reads with {args:?}");
    assert_eq!(dis, fastq(discarded), "discarded reads with {args:?}");
}

#[test]
fn discards_only_empty_reads_by_default() {
    check(&[], &[R1, R2, R3, R4, R5], &[E]);
}

#[test]
fn filters_on_expected_errors() {
    // r1 has 0.001 expected errors, r5 0.0003
    check(&["--prob", "0.0004"], &[R5], &[R1, R2, R3, R4, E]);
}

#[test]
fn filters_on_expected_errors_per_base() {
    check(&["--prob-rate", "0.005"], &[R1, R5], &[R2, R3, R4, E]);
}

#[test]
fn filters_on_length() {
    check(&["--min", "4", "--max", "9"], &[R4], &[R1, R2, R3, R5, E]);
}

#[test]
fn filters_on_ns_ignoring_case() {
    check(&["--bad", "1"], &[R1, R2, R4, R5], &[R3, E]);
}

#[test]
fn strips_both_ends() {
    // discarded reads are written trimmed
    check(
        &["--strip-left", "2", "--strip-right", "4"],
        &[
            ("r1", "GTAC", "IIII"),
            ("r2", "GTAC", "III+"),
            ("r3;size=3", "GTNN", "5555"),
        ],
        &[("r4;size=3 desc", "", ""), ("r5;ee=0.5;", "", ""), E],
    );
}

#[test]
fn truncates_to_length_discarding_shorter_reads() {
    check(
        &["--trunc-len", "5"],
        &[
            ("r1", "ACGTA", "IIIII"),
            ("r2", "ACGTA", "IIIII"),
            ("r3;size=3", "ACGTN", "55555"),
            R4,
        ],
        &[R5, E],
    );
}

#[test]
fn truncates_to_length_keeping_shorter_reads() {
    check(
        &["--trunc-keep", "5"],
        &[
            ("r1", "ACGTA", "IIIII"),
            ("r2", "ACGTA", "IIIII"),
            ("r3;size=3", "ACGTN", "55555"),
            R4,
            R5,
        ],
        &[E],
    );
}

#[test]
fn truncates_at_low_quality() {
    // truncated before the first base of quality 10 or less
    check(
        &["--trunc-qual", "10"],
        &[R1, ("r2", "ACGTA", "IIIII"), R3, R5],
        &[("r4;size=3 desc", "", ""), E],
    );
}

#[test]
fn filters_on_minimum_quality() {
    // r3's bases have quality 20, which is enough
    check(&["--min-qual", "20"], &[R1, R3, R5], &[R2, R4, E]);
}

#[test]
fn truncates_on_expected_errors() {
    check(
        &["--trunc-prob", "0.02"],
        &[R1, ("r2", "ACGTA", "IIIII"), ("r3;size=3", "AC", "55"), R5],
        &[("r4;size=3 desc", "", ""), E],
    );
}

#[test]
fn truncates_on_expected_errors_per_base() {
    check(
        &["--trunc-rate", "0.01"],
        &[
            R1,
            ("r2", "ACGTA", "IIIII"),
            ("r3;size=3", "ACGTN", "55555"),
            R5,
        ],
        &[("r4;size=3 desc", "", ""), E],
    );
}

#[test]
fn filters_on_abundance() {
    // a description after ";size=3" hides it, so r4 counts as size 1
    check(&["--min-size", "2"], &[R3], &[R1, R2, R4, R5, E]);
    check(&["--max-size", "2"], &[R1, R2, R4, R5], &[R3, E]);
}

#[test]
fn annotates_expected_errors() {
    // an existing ee= annotation is replaced
    let (out, _) = qualf(&["--ee"], INPUT);
    let expected = fastq(&[
        ("r1;ee=0.001000", R1.1, R1.2),
        ("r2;ee=0.5005", R2.1, R2.2),
        ("r3;size=3;ee=0.10000", R3.1, R3.2),
        ("r4;size=3 desc;ee=3.1548", R4.1, R4.2),
        ("r5;ee=0.0003000", R5.1, R5.2),
    ]);
    assert_eq!(out, expected);
}

#[test]
fn writes_fasta() {
    let (out, dis) = qualf(&["--fasta"], INPUT);
    let expected = ">r1\nACGTACGTAC\n>r2\nACGTACGTAC\n>r3;size=3\nACGTNNACGT\n\
                    >r4;size=3 desc\nacgtn\n>r5;ee=0.5;\nACG\n";
    assert_eq!(out, expected);
    assert_eq!(dis, ">e\n\n");
}

#[test]
fn reads_multiline_and_crlf_records() {
    let input = "@a\nACG\nTA\n+\nIII\nII\n@b\r\nAC\r\n+\r\nII\r\n";
    let (out, _) = qualf(&[], input);
    assert_eq!(out, fastq(&[("a", "ACGTA", "IIIII"), ("b", "AC", "II")]));
}

#[test]
fn drops_discarded_reads_without_discard_file() {
    let output = run(&["--min", "4"], INPUT);
    assert!(output.status.success());
    let out = String::from_utf8(output.stdout).unwrap();
    assert_eq!(out, fastq(&[R1, R2, R3, R4]));
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, INPUT).unwrap();
    let result = run(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ],
        "",
    );
    assert!(result.status.success());
    assert!(result.stdout.is_empty());
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(written, fastq(&[R1, R2, R3, R4, R5]));
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--prob", "0"][..],
        &["--prob-rate", "-1"],
        &["--trunc-qual", "94"],
        &["--min", "0"],
        &["--trunc-len", "0"],
        &["--min-size", "0"],
        &["--unknown"],
    ] {
        assert!(!run(args, INPUT).status.success(), "{args:?} should fail");
    }
}

#[test]
fn rejects_malformed_input() {
    for input in ["r1\nACGT\n+\nIIII\n", "@r1\nACGT\n+\nIII\n", "@r1\nACGT\n"] {
        assert!(!run(&[], input).status.success(), "{input:?} should fail");
    }
}
