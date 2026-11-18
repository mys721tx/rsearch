//! End-to-end tests of the derep binary.
//!
//! Unless noted otherwise, expected outputs are those of
//! `vsearch --fastx_uniques --sizein --sizeout --minseqlength 1` on the same
//! input, so that derep keeps vsearch's dereplication behaviour.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs derep with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_derep"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("derep should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs derep and returns its standard output, which must succeed.
fn derep(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "derep failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-derep-{}-{name}", std::process::id()))
}

#[test]
fn merges_identical_sequences() {
    let input = fasta(&[("a", "ACGT"), ("b", "ACGT"), ("c", "TTTT"), ("d", "acgt")]);
    let expected = fasta(&[("a;size=3", "ACGT"), ("c;size=1", "TTTT")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn ignores_case_and_u_versus_t() {
    // the first sequence is written as it appears in the input
    let input = fasta(&[("a", "ACGU"), ("b", "acgt")]);
    assert_eq!(derep(&[], &input), fasta(&[("a;size=2", "ACGU")]));
}

#[test]
fn keeps_ambiguity_codes_distinct() {
    let input = fasta(&[("a", "ACGTN"), ("b", "ACGTR"), ("c", "ACGTn")]);
    let expected = fasta(&[("a;size=2", "ACGTN"), ("b;size=1", "ACGTR")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn adds_up_input_sizes() {
    let input = fasta(&[("a;size=2", "AAAA"), ("b;size=3", "CCCC"), ("c", "AAAA")]);
    let expected = fasta(&[("a;size=3", "AAAA"), ("b;size=3", "CCCC")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn sorts_by_decreasing_size_then_label() {
    let input = fasta(&[
        ("b", "CCCC"),
        ("a", "GGGG"),
        ("c", "AAAA"),
        ("d", "TTTT"),
        ("e", "TTTT"),
    ]);
    let expected = fasta(&[
        ("d;size=2", "TTTT"),
        ("a;size=1", "GGGG"),
        ("b;size=1", "CCCC"),
        ("c;size=1", "AAAA"),
    ]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn compares_whole_labels_to_break_ties() {
    let input = fasta(&[("b;x=1", "AAAA"), ("b", "CCCC"), ("b;a=1", "GGGG")]);
    let expected = fasta(&[
        ("b;size=1", "CCCC"),
        ("b;a=1;size=1", "GGGG"),
        ("b;x=1;size=1", "AAAA"),
    ]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn breaks_ties_on_label_by_first_appearance() {
    let input = fasta(&[("x", "CCCC"), ("x", "AAAA")]);
    let expected = fasta(&[("x;size=1", "CCCC"), ("x;size=1", "AAAA")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn unsorted_keeps_order_of_first_appearance() {
    // vsearch --derep_smallmem
    let input = fasta(&[("d", "TTTT"), ("c", "GGGG"), ("a", "TTTT")]);
    let expected = fasta(&[("d;size=2", "TTTT"), ("c;size=1", "GGGG")]);
    assert_eq!(derep(&["--unsorted"], &input), expected);
}

#[test]
fn merges_reverse_complements_only_with_strand() {
    let input = fasta(&[("a", "AACC"), ("b", "GGTT"), ("c", "GGTT")]);
    let plus = fasta(&[("b;size=2", "GGTT"), ("a;size=1", "AACC")]);
    let both = fasta(&[("a;size=3", "AACC")]);
    assert_eq!(derep(&[], &input), plus);
    assert_eq!(derep(&["--strand"], &input), both);
    assert_eq!(derep(&["-s"], &input), both);
}

#[test]
fn prefers_the_plus_strand_match() {
    // GGTT is both a copy of the second sequence and the reverse
    // complement of the first; the identical sequence wins
    let input = fasta(&[("a", "AACC"), ("b", "GGTT"), ("c", "GGTT")]);
    let output = derep(&["--unsorted", "--strand"], &input);
    assert_eq!(output, fasta(&[("a;size=3", "AACC")]));
    let input = fasta(&[("b", "GGTT"), ("a", "AACC"), ("c", "GGTT")]);
    let output = derep(&["--unsorted", "--strand"], &input);
    assert_eq!(output, fasta(&[("b;size=3", "GGTT")]));
}

#[test]
fn filters_by_size_and_relabels_written_sequences() {
    let input = fasta(&[
        ("a;size=5", "AAAA"),
        ("b;size=3", "CCCC"),
        ("c;size=2", "GGGG"),
        ("d", "TTTT"),
    ]);
    let expected = fasta(&[("U1;size=3", "CCCC"), ("U2;size=2", "GGGG")]);
    let args = ["--min", "2", "--max", "4", "--relabel", "U"];
    assert_eq!(derep(&args, &input), expected);
}

#[test]
fn writes_only_the_top_sequences() {
    let input = fasta(&[
        ("a;size=5", "AAAA"),
        ("b;size=3", "CCCC"),
        ("c;size=2", "GGGG"),
    ]);
    let expected = fasta(&[("a;size=5", "AAAA"), ("b;size=3", "CCCC")]);
    assert_eq!(derep(&["--top", "2"], &input), expected);
}

#[test]
fn drops_empty_sequences() {
    let input = fasta(&[("a", "AAAA"), ("e", ""), ("b;size=3", "CCCC")]);
    let expected = fasta(&[("b;size=3", "CCCC"), ("a;size=1", "AAAA")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn joins_multiline_sequences() {
    let input = ">a\nAC\nGT\n>b\nACGT\n";
    assert_eq!(derep(&[], input), fasta(&[("a;size=2", "ACGT")]));
}

#[test]
fn replaces_only_the_size_annotation() {
    // the description is dropped, other annotations are kept, and the
    // first valid size annotation is replaced by the new size at the end
    for (label, expected) in [
        ("sample=1;foo;size=3 desc", "sample=1;foo;size=4"),
        ("a;size=3;x=1", "a;x=1;size=4"),
        ("a;size=2;size=5", "a;size=5;size=3"),
        ("size=3", ";size=4"),
        ("a;size=abc", "a;size=abc;size=2"),
        ("asize=2", "asize=2;size=2"),
    ] {
        let input = fasta(&[(label, "ACGT"), ("z", "ACGT")]);
        assert_eq!(derep(&[], &input), fasta(&[(expected, "ACGT")]), "{label}");
    }
}

#[test]
fn rejects_zero_sizes() {
    let output = run(&[], &fasta(&[("a;size=0", "ACGT")]));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
}

#[test]
fn wraps_fasta_lines() {
    let seq = "ACGT".repeat(60);
    let input = fasta(&[("a", &seq), ("b", &seq)]);
    // vsearch's default width of 80
    let wrapped = format!(
        ">a;size=2\n{}\n{}\n{}\n",
        &seq[..80],
        &seq[80..160],
        &seq[160..]
    );
    assert_eq!(derep(&["--fasta-width", "80"], &input), wrapped);
    // vsearch --fasta_width 0
    assert_eq!(derep(&[], &input), fasta(&[("a;size=2", &seq)]));
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, fasta(&[("a", "ACGT"), ("b", "ACGT")])).unwrap();
    let stdout = derep(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ],
        "",
    );
    assert_eq!(stdout, "");
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(written, fasta(&[("a;size=2", "ACGT")]));
}

#[test]
fn maps_members_to_written_sequences() {
    let map = temp_path("map.tsv");
    let input = fasta(&[
        ("a;size=2", "ACGT"),
        ("b", "TTTT"),
        ("c", "acgt"),
        ("d", "GGGG"),
    ]);
    let output = derep(
        &[
            "--min",
            "2",
            "--relabel",
            "U",
            "--map",
            map.to_str().unwrap(),
        ],
        &input,
    );
    let table = std::fs::read_to_string(&map).unwrap();
    std::fs::remove_file(map).unwrap();
    assert_eq!(output, fasta(&[("U1;size=3", "ACGT")]));
    // derep's own format: unlike vsearch --tabbedout, members of sequences
    // that are not written are left out, and labels drop their sizes
    assert_eq!(table, "a\tU1\nc\tU1\n");
}

#[test]
fn merges_fastq_qualities() {
    let input = fastq(&[("a", "ACGT", "I#5+"), ("b", "ACGT", "#I5,")]);
    assert_eq!(derep(&[], &input), fastq(&[("a;size=2", "ACGT", "&&5+")]));
}

#[test]
fn weights_fastq_qualities_by_size() {
    let input = fastq(&[
        ("a;size=3", "ACGT", "5555"),
        ("b;size=2", "acgu", "++++"),
        ("c", "ACGT", "????"),
    ]);
    assert_eq!(derep(&[], &input), fastq(&[("a;size=6", "ACGT", "....")]));
}

#[test]
fn reverses_qualities_of_reverse_strand_reads() {
    // vsearch --fastx_uniques --strand both
    let input = fastq(&[("a", "AACC", "ABCD"), ("b", "GGTT", "IIII")]);
    let expected = fastq(&[("a;size=2", "AACC", "CDEE")]);
    assert_eq!(derep(&["--strand"], &input), expected);
}

#[test]
fn drops_empty_fastq_reads() {
    let input = fastq(&[("a", "ACGT", "IIII"), ("e", "", ""), ("b", "ACGT", "IIII")]);
    assert_eq!(derep(&[], &input), fastq(&[("a;size=2", "ACGT", "IIII")]));
}

#[test]
fn splits_fasta_labels_at_any_whitespace() {
    // the label is the id, which ends at a tab as well as a space
    let input = ">a\tdesc\nACGT\n>b desc\nACGT\n";
    assert_eq!(derep(&[], input), fasta(&[("a;size=2", "ACGT")]));
}

#[test]
fn raises_quality_zero_to_one() {
    let input = fastq(&[("a", "ACGT", "!!+I")]);
    assert_eq!(derep(&[], &input), fastq(&[("a;size=1", "ACGT", "\"\"+I")]));
}

#[test]
fn rejects_invalid_arguments() {
    for args in [&["--min", "x"][..], &["--top", "-1"], &["--unknown"]] {
        assert!(!run(args, "").status.success(), "{args:?} should fail");
    }
}

#[test]
fn accepts_empty_input() {
    assert_eq!(derep(&[], ""), "");
}
