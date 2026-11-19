//! End-to-end tests of the subsmp binary.
//!
//! Expected outputs are those of `vsearch --fastx_subsample --notrunclabels
//! --minseqlength 0 --fasta_width 0` with `--randseed` (`--seed`) and
//! `--sample_size` (`--count`) or `--sample_pct` (`--pct`) on the same
//! input. The reads drawn are written in input order, and the others with
//! `--fastaout_discarded` (`--discarded`). `--by-size` is vsearch's
//! `--sizein --sizeout`. Each record counts as its abundance, and is
//! written with the number of its reads drawn. Without a seed the draw is
//! checked only for its size.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs subsmp with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_subsmp"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("subsmp should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs subsmp and returns its standard output, which must succeed.
fn subsmp(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "subsmp failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-subsmp-{}-{name}", std::process::id()))
}

/// Six records with abundances 5, 1, 1 (one hidden by a description), 2,
/// 10 and 1 (an empty sequence): six records, or twenty reads by size.
fn records() -> String {
    fasta(&[
        ("a;size=5", "A"),
        ("b", "C"),
        ("c;size=3 x", "G"),
        ("d;size=2", "T"),
        ("e;size=10", "AC"),
        ("f", ""),
    ])
}

#[test]
fn draws_a_number_of_records() {
    // headers are written unchanged
    let expected = fasta(&[("a;size=5", "A"), ("b", "C"), ("d;size=2", "T")]);
    assert_eq!(
        subsmp(&["--seed", "1", "--count", "3"], &records()),
        expected
    );
    let expected = fasta(&[("a;size=5", "A"), ("c;size=3 x", "G")]);
    assert_eq!(
        subsmp(&["--seed", "-1", "--count", "2"], &records()),
        expected
    );
}

#[test]
fn draws_a_share_of_records() {
    // the number drawn is rounded down: 50% of 6 is 3, 33.3% is 1
    let expected = fasta(&[("a;size=5", "A"), ("b", "C"), ("d;size=2", "T")]);
    assert_eq!(
        subsmp(&["--seed", "1", "--pct", "50"], &records()),
        expected
    );
    let expected = fasta(&[("a;size=5", "A")]);
    assert_eq!(
        subsmp(&["--seed", "1", "--pct", "33.3"], &records()),
        expected
    );
}

#[test]
fn draws_reads_by_abundance() {
    for (args, expected) in [
        (
            &["--seed", "1", "--count", "3"][..],
            &[("a;size=2", "A"), ("d;size=1", "T")][..],
        ),
        (
            &["--seed", "42", "--count", "10"],
            &[
                ("a;size=1", "A"),
                ("b;size=1", "C"),
                ("d;size=2", "T"),
                ("e;size=5", "AC"),
                ("f;size=1", ""),
            ],
        ),
        (
            &["--seed", "1", "--pct", "50"],
            &[("a;size=4", "A"), ("d;size=1", "T"), ("e;size=5", "AC")],
        ),
        (
            // 33.3% of 20 reads is 6
            &["--seed", "7", "--pct", "33.3"],
            &[
                ("a;size=2", "A"),
                ("b;size=1", "C"),
                ("e;size=2", "AC"),
                ("f;size=1", ""),
            ],
        ),
    ] {
        let mut all = args.to_vec();
        all.push("--by-size");
        assert_eq!(subsmp(&all, &records()), fasta(expected), "{args:?}");
    }
}

#[test]
fn writes_the_reads_not_drawn() {
    let path = temp_path("discarded.fa");
    let args = [
        "--seed",
        "1",
        "--count",
        "3",
        "--by-size",
        "--discarded",
        path.to_str().unwrap(),
    ];
    let kept = subsmp(&args, &records());
    let discarded = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert_eq!(kept, fasta(&[("a;size=2", "A"), ("d;size=1", "T")]));
    // the rest of each abundance
    let expected = fasta(&[
        ("a;size=3", "A"),
        ("b;size=1", "C"),
        ("c;size=3 x;size=1", "G"),
        ("d;size=1", "T"),
        ("e;size=10", "AC"),
        ("f;size=1", ""),
    ]);
    assert_eq!(discarded, expected);
}

#[test]
fn draws_fastq() {
    let path = temp_path("discarded.fq");
    let input = fastq(&[
        ("q;size=2", "AC", "II"),
        ("r", "G", "5"),
        ("s;size=3", "TT", "#I"),
    ]);
    let args = [
        "--seed",
        "3",
        "--count",
        "3",
        "--by-size",
        "--discarded",
        path.to_str().unwrap(),
    ];
    let kept = subsmp(&args, &input);
    let discarded = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert_eq!(
        kept,
        fastq(&[("q;size=1", "AC", "II"), ("s;size=2", "TT", "#I")])
    );
    assert_eq!(
        discarded,
        fastq(&[
            ("q;size=1", "AC", "II"),
            ("r;size=1", "G", "5"),
            ("s;size=1", "TT", "#I"),
        ])
    );
}

#[test]
fn draws_everything_when_asked_for_all() {
    assert_eq!(
        subsmp(&["--seed", "1", "--count", "6"], &records()),
        records()
    );
    assert_eq!(
        subsmp(&["--seed", "1", "--pct", "100"], &records()),
        records()
    );
}

#[test]
fn refuses_more_than_there_is_unless_allowed() {
    let output = run(&["--seed", "1", "--count", "7"], &records());
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Cannot subsample 7 reads from a sample of 6"),
        "{stderr}"
    );
    let args = ["--seed", "1", "--count", "7", "--allow-fewer"];
    assert_eq!(subsmp(&args, &records()), records());
}

#[test]
fn draws_from_empty_input() {
    assert!(!run(&["--count", "1"], "").status.success());
    assert_eq!(subsmp(&["--count", "1", "--allow-fewer"], ""), "");
    assert_eq!(subsmp(&["--pct", "50"], ""), "");
}

#[test]
fn draws_with_a_system_seed_without_one() {
    for args in [&["--count", "3"][..], &["--count", "3", "--seed", "0"]] {
        let output = subsmp(args, &records());
        assert_eq!(output.matches('>').count(), 3, "{args:?}");
    }
}

#[test]
fn needs_one_size_of_draw() {
    // vsearch: specify either --sample_pct or --sample_size, a size of 0
    // counting as none
    for args in [
        &[][..],
        &["--count", "2", "--pct", "50"],
        &["--count", "0"],
        &["--count", "-1"],
        &["--pct", "0"],
        &["--pct", "-1"],
        &["--pct", "100.1"],
        &["--count", "2", "--seed", "x"],
    ] {
        assert!(
            !run(args, &records()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn rejects_invalid_abundances_before_writing() {
    // vsearch reads every abundance, with or without --sizein
    let input = fasta(&[("b", "C"), ("a;size=0", "A")]);
    for args in [&["--count", "1"][..], &["--count", "1", "--by-size"]] {
        let output = run(args, &input);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
    }
}

#[test]
fn rejects_invalid_input_before_writing() {
    let output = run(&["--count", "1"], ">a\nACGT\n>b\nAC-T\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, records()).unwrap();
    let stdout = subsmp(
        &[
            "--seed",
            "1",
            "--count",
            "3",
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
    assert_eq!(
        written,
        fasta(&[("a;size=5", "A"), ("b", "C"), ("d;size=2", "T")])
    );
}
