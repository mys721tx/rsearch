//! End-to-end tests of the subseq binary.
//!
//! Expected outputs are those of `vsearch --fastx_getsubseq
//! --notrunclabels` with `--subseq_start` (`--start`) and `--subseq_end`
//! (`--end`), selecting every record (`--label_substr_match --label ''`).
//! Each read's region from start to end is written, both counted from 1 and
//! included. A region past the end of a read leaves it empty. Selecting
//! records is getseq's job, and `--fastx_getsubseq --label L` is `getseq
//! --label L | subseq`.

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

/// Runs subseq with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_subseq"), args, input)
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

/// Runs subseq and returns its standard output, which must succeed.
fn subseq(args: &[&str], input: &str) -> String {
    stdout_of("subseq", run(args, input))
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
    std::env::temp_dir().join(format!("rsearch-subseq-{}-{name}", std::process::id()))
}

/// Reads of 8, 2 and 5 bases.
fn reads() -> String {
    fastq(&[
        ("a desc", "ACGTACGT", "!\"#$%&'("),
        ("b;size=2", "AC", "II"),
        ("c", "ACGTA", "55555"),
    ])
}

/// The three reads cut to the given sequences and scores.
fn cut(regions: [(&str, &str); 3]) -> String {
    let labels = ["a desc", "b;size=2", "c"];
    let records: Vec<_> = (0..3)
        .map(|i| (labels[i], regions[i].0, regions[i].1))
        .collect();
    fastq(&records)
}

#[test]
fn keeps_a_region() {
    let expected = cut([("GTAC", "#$%&"), ("", ""), ("GTA", "555")]);
    assert_eq!(subseq(&["--start", "3", "--end", "6"], &reads()), expected);
}

#[test]
fn keeps_one_base() {
    let expected = cut([("G", "#"), ("", ""), ("G", "5")]);
    assert_eq!(subseq(&["--start", "3", "--end", "3"], &reads()), expected);
}

#[test]
fn keeps_to_the_end_by_default() {
    let expected = cut([("GTACGT", "#$%&'("), ("", ""), ("GTA", "555")]);
    assert_eq!(subseq(&["--start", "3"], &reads()), expected);
}

#[test]
fn keeps_from_the_start_by_default() {
    let expected = cut([("AC", "!\""), ("AC", "II"), ("AC", "55")]);
    assert_eq!(subseq(&["--end", "2"], &reads()), expected);
    assert_eq!(subseq(&[], &reads()), reads());
}

#[test]
fn empties_reads_shorter_than_the_start() {
    let expected = cut([("CGT", "&'("), ("", ""), ("", "")]);
    assert_eq!(
        subseq(&["--start", "6", "--end", "100"], &reads()),
        expected
    );
    let expected = cut([("", ""), ("", ""), ("", "")]);
    assert_eq!(subseq(&["--start", "9", "--end", "9"], &reads()), expected);
}

#[test]
fn cuts_fasta() {
    // the sequence as read, lines joined and invalid characters stripped
    let output = run(
        &["--start", "2", "--end", "6"],
        ">a x\nACGT\nACGT\n>b\nAC1GT\n",
    );
    assert!(output.status.success());
    let expected = fasta(&[("a x", "CGTAC"), ("b", "CGT")]);
    assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn rejects_invalid_regions() {
    // vsearch: at least 1, and the start no greater than the end
    for args in [
        &["--start", "0", "--end", "3"][..],
        &["--start", "-1"],
        &["--end", "0"],
        &["--start", "4", "--end", "3"],
        &["--start", "x"],
    ] {
        assert!(
            !run(args, &reads()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn rejects_invalid_abundances() {
    // vsearch reads the abundance of every record it writes
    let output = run(&["--start", "2"], &fasta(&[("a;size=0", "ACGT")]));
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
    let input = fasta(&[("a;size=1", "ACGT"), ("z;size=9223372036854775808", "AC")]);
    let output = run(&["--start", "2"], &input);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a;size=1", "CGT")])
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("range"));
}

#[test]
fn rejects_invalid_input() {
    let output = run(&["--start", "2"], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character 'I'"),
        "{stderr}"
    );
}

#[test]
fn accepts_empty_input() {
    assert_eq!(subseq(&["--start", "2"], ""), "");
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, reads()).unwrap();
    let stdout = subseq(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--end",
            "2",
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(written, cut([("AC", "!\""), ("AC", "II"), ("AC", "55")]));
}

#[test]
fn pipeline_fastx_getsubseq() {
    // vsearch --fastx_getsubseq in --notrunclabels --label 'a desc'
    // --subseq_start 3 --subseq_end 6 --fastqout out
    let getseq = env!("CARGO_BIN_EXE_getseq");
    let selected = stdout_of("getseq", run_bin(getseq, &["--label", "a desc"], &reads()));
    let expected = fastq(&[("a desc", "GTAC", "#$%&")]);
    assert_eq!(subseq(&["--start", "3", "--end", "6"], &selected), expected);
}
