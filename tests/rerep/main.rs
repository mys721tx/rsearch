//! End-to-end tests of the rerep binary.
//!
//! Expected outputs are those of `vsearch --rereplicate --notrunclabels
//! --fasta_width 0 --sizeout` on the same input. Each read is written as
//! many times as its abundance (`;size=`), each copy with an abundance of
//! 1. `--no-sizeout` is vsearch without `--sizeout`, which removes the
//! abundance. vsearch reads only FASTA; FASTQ copies keep their scores, and
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
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs rerep with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_rerep"), args, input)
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

/// Runs rerep and returns its standard output, which must succeed.
fn rerep(args: &[&str], input: &str) -> String {
    stdout_of("rerep", run(args, input))
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
    std::env::temp_dir().join(format!("rsearch-rerep-{}-{name}", std::process::id()))
}

/// Reads with and without abundances: one of 3, one without, one whose
/// abundance sits among other annotations and a description, one of 1, an
/// empty one of 2, and one whose abundance a description hides.
fn reads() -> String {
    fasta(&[
        ("a;size=3", "ACGT"),
        ("b", "GG"),
        ("c;size=2;x=1 desc", "TT"),
        ("d;size=1", "CC"),
        ("e;size=2", ""),
        ("f;size=2 x", "AA"),
    ])
}

#[test]
fn writes_each_read_as_many_times_as_its_abundance() {
    // a read without an abundance is written once, and the abundance of
    // every copy is 1
    let expected = fasta(&[
        ("a;size=1", "ACGT"),
        ("a;size=1", "ACGT"),
        ("a;size=1", "ACGT"),
        ("b;size=1", "GG"),
        ("c;x=1 desc;size=1", "TT"),
        ("c;x=1 desc;size=1", "TT"),
        ("d;size=1", "CC"),
        ("e;size=1", ""),
        ("e;size=1", ""),
        ("f;size=2 x;size=1", "AA"),
    ]);
    assert_eq!(rerep(&[], &reads()), expected);
}

#[test]
fn removes_abundances_without_sizeout() {
    // vsearch --rereplicate without --sizeout
    let expected = fasta(&[
        ("a", "ACGT"),
        ("a", "ACGT"),
        ("a", "ACGT"),
        ("b", "GG"),
        ("c;x=1 desc", "TT"),
        ("c;x=1 desc", "TT"),
        ("d", "CC"),
        ("e", ""),
        ("e", ""),
        ("f;size=2 x", "AA"),
    ]);
    assert_eq!(rerep(&["--no-sizeout"], &reads()), expected);
}

#[test]
fn warns_about_missing_abundances() {
    let output = run(&[], &reads());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Missing abundance information for some input sequences, assumed 1"),
        "{stderr}"
    );
    let output = run(&[], &fasta(&[("a;size=2", "AC")]));
    assert!(output.stderr.is_empty());
}

#[test]
fn keeps_whole_headers() {
    // vsearch --notrunclabels
    let input = fasta(&[("a\tt;size=2", "AC")]);
    let expected = fasta(&[("a\tt;size=1", "AC"), ("a\tt;size=1", "AC")]);
    assert_eq!(rerep(&[], &input), expected);
}

#[test]
fn unwraps_fasta_and_strips_invalid_characters() {
    let output = run(&[], ">w;size=2\nAC1GT\nAC\n");
    assert!(output.status.success());
    let expected = fasta(&[("w;size=1", "ACGTAC"), ("w;size=1", "ACGTAC")]);
    assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn rejects_zero_abundance() {
    let output = run(&[], &fasta(&[("z;size=0", "AC")]));
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
}

#[test]
fn rejects_abundance_out_of_range() {
    // the reads before are written
    let input = fasta(&[("a;size=1", "AC"), ("z;size=9223372036854775808", "AC")]);
    let output = run(&[], &input);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a;size=1", "AC")])
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("range"));
}

#[test]
fn copies_fastq_scores() {
    let input = fastq(&[("q;size=2", "AC", "I#"), ("r", "G", "5")]);
    let expected = fastq(&[
        ("q;size=1", "AC", "I#"),
        ("q;size=1", "AC", "I#"),
        ("r;size=1", "G", "5"),
    ]);
    assert_eq!(rerep(&[], &input), expected);
    let expected = fastq(&[("q", "AC", "I#"), ("q", "AC", "I#"), ("r", "G", "5")]);
    assert_eq!(rerep(&["--no-sizeout"], &input), expected);
}

#[test]
fn matches_vsearch_on_fastq_written_as_fasta() {
    // rerep | conv --fasta --fasta-width 0 is vsearch --rereplicate on the
    // reads as FASTA
    let input = fastq(&[
        ("a;size=3", "ACGT", "IIII"),
        ("b", "GG", "##"),
        ("c;size=2;x=1 desc", "TT", "5+"),
    ]);
    let copies = rerep(&[], &input);
    let conv = env!("CARGO_BIN_EXE_conv");
    let output = stdout_of(
        "conv",
        run_bin(conv, &["--fasta", "--fasta-width", "0"], &copies),
    );
    let expected = fasta(&[
        ("a;size=1", "ACGT"),
        ("a;size=1", "ACGT"),
        ("a;size=1", "ACGT"),
        ("b;size=1", "GG"),
        ("c;x=1 desc;size=1", "TT"),
        ("c;x=1 desc;size=1", "TT"),
    ]);
    assert_eq!(output, expected);
}

#[test]
fn accepts_empty_input() {
    assert_eq!(rerep(&[], ""), "");
}

#[test]
fn rejects_invalid_input() {
    let output = run(&[], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character 'I'"),
        "{stderr}"
    );
}

#[test]
fn rejects_invalid_arguments() {
    for args in [&["--unknown"][..], &["--sizeout"], &["--relabel", "R"]] {
        assert!(
            !run(args, &reads()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, fasta(&[("a;size=2", "AC")])).unwrap();
    let stdout = rerep(
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
    assert_eq!(written, fasta(&[("a;size=1", "AC"), ("a;size=1", "AC")]));
}

#[test]
fn pipeline_dereplication_and_rereplication() {
    // vsearch --fastx_uniques in --sizein --sizeout --notrunclabels
    // --minseqlength 0, then --rereplicate --sizeout
    let input = fasta(&[
        ("a", "ACGT"),
        ("b", "TTTT"),
        ("c", "acgt"),
        ("d;size=2", "GG"),
        ("e", "TTTT"),
    ]);
    let derep = env!("CARGO_BIN_EXE_derep");
    let dereplicated = stdout_of("derep", run_bin(derep, &[], &input));
    let expected = fasta(&[
        ("a;size=1", "ACGT"),
        ("a;size=1", "ACGT"),
        ("b;size=1", "TTTT"),
        ("b;size=1", "TTTT"),
        ("d;size=1", "GG"),
        ("d;size=1", "GG"),
    ]);
    assert_eq!(rerep(&[], &dereplicated), expected);
}
