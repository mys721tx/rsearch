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
    // a file name unique to the arguments
    let name: String = args
        .join("_")
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    let discard = temp_path(&format!("discard-{name}"));
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

/// FASTA reads: an abundance, two Ns with a description, stale size and
/// length annotations, and an empty sequence.
const FASTA: &str =
    ">f1;size=3\nACGTACGTAC\n>f2 desc\nACGNNACGT\n>f3;size=2;length=5\nACG\n>f4\n\n";

#[test]
fn filters_fasta_by_length_ns_and_size() {
    // FASTA output, with the length, N and abundance filters
    let (out, dis) = qualf(&["--min", "4", "--bad", "1"], FASTA);
    assert_eq!(out, ">f1;size=3\nACGTACGTAC\n");
    assert_eq!(
        dis,
        ">f2 desc\nACGNNACGT\n>f3;size=2;length=5\nACG\n>f4\n\n"
    );
    let (out, _) = qualf(&["--min-size", "2"], FASTA);
    assert_eq!(out, ">f1;size=3\nACGTACGTAC\n>f3;size=2;length=5\nACG\n");
}

#[test]
fn rejects_malformed_input() {
    for input in ["r1\nACGT\n+\nIIII\n", "@r1\nACGT\n+\nIII\n", "@r1\nACGT\n"] {
        assert!(!run(&[], input).status.success(), "{input:?} should fail");
    }
}

#[test]
fn rejects_quality_options_with_fasta() {
    for option in [
        &["--prob", "1"][..],
        &["--prob-rate", "0.1"],
        &["--min-qual", "5"],
        &["--trunc-qual", "5"],
        &["--trunc-prob", "1"],
        &["--trunc-rate", "0.1"],
        &["--ee"],
        &["--ascii", "64"],
        &["--qmin", "1"],
        &["--qmax", "40"],
    ] {
        let output = run(option, FASTA);
        assert!(!output.status.success(), "{option:?} should fail");
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn annotates_size_and_length() {
    // the old annotations are replaced
    let (out, _) = qualf(&["--sizeout", "--lengthout"], FASTA);
    let expected = ">f1;size=3;length=10\nACGTACGTAC\n>f2 desc;size=1;length=9\nACGNNACGT\n\
                    >f3;size=2;length=3\nACG\n";
    assert_eq!(out, expected);
}

#[test]
fn removes_annotations_and_adds_sample_and_suffix() {
    let (out, _) = qualf(
        &[
            "--xsize",
            "--xlength",
            "--sample",
            "S",
            "--label-suffix",
            "_x",
        ],
        FASTA,
    );
    let expected =
        ">f1_x;sample=S\nACGTACGTAC\n>f2 desc_x;sample=S\nACGNNACGT\n>f3_x;sample=S\nACG\n";
    assert_eq!(out, expected);

    let input = "@q;ee=9;size=4\nACGT\n+\nIIII\n";
    let (out, _) = qualf(&["--xee", "--sizeout"], input);
    assert_eq!(out, fastq(&[("q;size=4", "ACGT", "IIII")]));
}

#[test]
fn relabels_by_position() {
    // kept and discarded sequences are numbered separately
    let (out, dis) = qualf(&["--relabel", "R", "--relabel-keep"], FASTA);
    let expected = ">R1 f1;size=3\nACGTACGTAC\n>R2 f2 desc\nACGNNACGT\n\
                    >R3 f3;size=2;length=5\nACG\n";
    assert_eq!(out, expected);
    assert_eq!(dis, ">R1 f4\n\n");
}

#[test]
fn relabels_by_sequence() {
    let (out, _) = qualf(&["--relabel-self", "--sizeout"], FASTA);
    let expected =
        ">ACGTACGTAC;size=3\nACGTACGTAC\n>ACGNNACGT;size=1\nACGNNACGT\n>ACG;size=2\nACG\n";
    assert_eq!(out, expected);

    let (out, _) = qualf(&["--relabel-sha1"], FASTA);
    let expected = ">e8e5136aebc65f03c77c95fb40e0e2502de7e9a8\nACGTACGTAC\n\
                    >d6a1ef23383fc7561ed765f9a6545c1f5b748b7e\nACGNNACGT\n\
                    >b88cad96aeb0ad3b4a10da7d0372a759067df35e\nACG\n";
    assert_eq!(out, expected);

    let (out, _) = qualf(&["--relabel-md5"], FASTA);
    let expected = ">45aff2fecf7615d56bc0567dffab9fa8\nACGTACGTAC\n\
                    >0a659c2752840dd0e3cd2aab30accdbf\nACGNNACGT\n\
                    >33f786e15eb427ffd3edec16cfdc0cd2\nACG\n";
    assert_eq!(out, expected);
}

#[test]
fn rejects_more_than_one_relabelling() {
    let output = run(&["--relabel", "R", "--relabel-md5"], FASTA);
    assert!(!output.status.success());
}

/// Paired reads: p1's reverse read fails --prob 0.5, p2's fails badly,
/// p3 passes in both directions.
const FORWARD: &str =
    "@p1\nACGTACGTAC\n+\nIIIIIIIIII\n@p2\nACGTACGTAC\n+\nIIIIIIIIII\n@p3\nACGT\n+\nIIII\n";
const REVERSE: &str =
    "@p1\nGTACGTACGT\n+\nIIIII+++++\n@p2\nGTACGTACGT\n+\n##########\n@p3\nACGT\n+\nIIII\n";

#[test]
fn filters_pairs_together() {
    let paths = ["rev.fq", "rev-out.fq", "rev-dis.fq"].map(temp_path);
    let [rev, rev_out, rev_dis] = paths.each_ref().map(|p| p.to_str().unwrap());
    std::fs::write(rev, REVERSE).unwrap();
    let args = [
        "--prob",
        "0.5",
        "--relabel",
        "P",
        "--reverse",
        rev,
        "--rev-out",
        rev_out,
        "--rev-discard",
        rev_dis,
    ];
    let (out, dis) = qualf(&args, FORWARD);
    let (kept_rev, dis_rev) = (
        std::fs::read_to_string(rev_out).unwrap(),
        std::fs::read_to_string(rev_dis).unwrap(),
    );
    for path in &paths {
        std::fs::remove_file(path).unwrap();
    }
    // a pair is discarded when either read is; both reads get the same number
    assert_eq!(out, fastq(&[("P1", "ACGT", "IIII")]));
    assert_eq!(kept_rev, fastq(&[("P1", "ACGT", "IIII")]));
    assert_eq!(
        dis,
        fastq(&[
            ("P1", "ACGTACGTAC", "IIIIIIIIII"),
            ("P2", "ACGTACGTAC", "IIIIIIIIII")
        ])
    );
    assert_eq!(
        dis_rev,
        fastq(&[
            ("P1", "GTACGTACGT", "IIIII+++++"),
            ("P2", "GTACGTACGT", "##########")
        ])
    );
}

#[test]
fn rejects_unpaired_reads() {
    let rev = temp_path("short-rev.fq");
    let rev_str = rev.to_str().unwrap();
    std::fs::write(&rev, "@p1\nACGT\n+\nIIII\n").unwrap();
    let more_forward = run(&["--reverse", rev_str], FORWARD);
    std::fs::write(&rev, format!("{REVERSE}@p4\nA\n+\nI\n")).unwrap();
    let more_reverse = run(&["--reverse", rev_str], FORWARD);
    std::fs::write(&rev, ">p1\nACGT\n").unwrap();
    let other_format = run(&["--reverse", rev_str], FORWARD);
    std::fs::remove_file(rev).unwrap();
    assert!(!more_forward.status.success());
    assert!(!more_reverse.status.success());
    assert!(!other_format.status.success());
}

#[test]
fn needs_reverse_reads_for_reverse_outputs() {
    assert!(!run(&["--rev-out", "/dev/null"], FORWARD).status.success());
    assert!(
        !run(&["--rev-discard", "/dev/null"], FORWARD)
            .status
            .success()
    );
}

/// Phred+64 reads: `a` is Q40 throughout, `b` starts with two scores of -5.
const PHRED64: &str = "@a\nACGT\n+\nhhhh\n@b\nACGT\n+\n;;hh\n";

#[test]
fn decodes_phred64() {
    // b's negative scores fall below the default --min-qual of 0
    let (out, dis) = qualf(&["--ascii", "64", "--qmin", "-5", "--ee"], PHRED64);
    assert_eq!(out, fastq(&[("a;ee=0.0004000", "ACGT", "hhhh")]));
    assert_eq!(dis, fastq(&[("b;ee=6.3248", "ACGT", ";;hh")]));
}

#[test]
fn rejects_scores_outside_the_range() {
    // with the default --qmin of 0, b's scores of -5 are an error
    let output = run(&["--ascii", "64"], PHRED64);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("below the minimum 0 in read b"), "{stderr}");

    // the last two bases are Q41
    let input = "@a\nACGTACGT\n+\nIIIIIIJJ\n";
    let output = run(&["--qmax", "40"], input);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("above the maximum 40 in read a"),
        "{stderr}"
    );

    // only the bases left after trimming are checked
    let (out, _) = qualf(&["--qmax", "40", "--trunc-keep", "6"], input);
    assert_eq!(out, fastq(&[("a", "ACGTAC", "IIIIII")]));
}

#[test]
fn rejects_invalid_quality_ranges() {
    for args in [
        &["--qmin", "10", "--qmax", "5"][..],
        &["--qmin", "-1"],
        &["--ascii", "64", "--qmax", "63"],
        &["--ascii", "50"],
    ] {
        assert!(!run(args, INPUT).status.success(), "{args:?} should fail");
    }
}

#[test]
fn wraps_fasta_lines() {
    let (out, _) = qualf(&["--fasta-width", "4"], ">w\nACGTACGTAC\n");
    assert_eq!(out, ">w\nACGT\nACGT\nAC\n");
    // FASTQ output is never wrapped
    let (out, _) = qualf(&["--fasta-width", "4"], INPUT);
    assert_eq!(out, fastq(&[R1, R2, R3, R4, R5]));
}

#[test]
fn prints_a_summary_unless_quiet() {
    let output = run(&["--trunc-keep", "6"], FORWARD);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "3 sequences kept (of which 2 truncated), 0 sequences discarded.\n"
    );
    let output = run(&["--min", "4"], FASTA);
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "2 sequences kept (of which 0 truncated), 2 sequences discarded.\n"
    );
    for quiet in ["-q", "--quiet"] {
        let output = run(&[quiet], FORWARD);
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
    }
}
