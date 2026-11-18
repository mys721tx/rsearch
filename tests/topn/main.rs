//! End-to-end tests of the topn binary.
//!
//! vsearch has no command of its own for this. `--topn N` keeps the first N
//! records a command writes. Expected outputs are therefore those of the
//! vsearch commands with `--topn`, compared with the rsearch pipeline that
//! ends with topn. Like vsearch, topn reads and checks its whole input.
//! Invalid input after the first N records thus still fails.

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

/// Runs topn with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_topn"), args, input)
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

/// Runs topn and returns its standard output, which must succeed.
fn topn(args: &[&str], input: &str) -> String {
    stdout_of("topn", run(args, input))
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
    std::env::temp_dir().join(format!("rsearch-topn-{}-{name}", std::process::id()))
}

/// Four records, one of them empty and one with a description.
fn records() -> String {
    fasta(&[
        ("a;size=3 desc", "ACGT"),
        ("b", ""),
        ("c", "GG"),
        ("d", "T"),
    ])
}

#[test]
fn keeps_the_first_records() {
    let expected = fasta(&[("a;size=3 desc", "ACGT"), ("b", "")]);
    assert_eq!(topn(&["2"], &records()), expected);
    assert_eq!(
        topn(&["1"], &records()),
        fasta(&[("a;size=3 desc", "ACGT")])
    );
}

#[test]
fn keeps_every_record_of_shorter_input() {
    assert_eq!(topn(&["4"], &records()), records());
    assert_eq!(topn(&["100"], &records()), records());
}

#[test]
fn keeps_fastq_records() {
    let input = fastq(&[("a", "ACGT", "II5I"), ("b", "TT", "#I"), ("c", "G", "5")]);
    let expected = fastq(&[("a", "ACGT", "II5I"), ("b", "TT", "#I")]);
    assert_eq!(topn(&["2"], &input), expected);
}

#[test]
fn unwraps_fasta() {
    let input = ">a\nAC\nGT\n>b\nTT\n";
    assert_eq!(topn(&["1"], input), fasta(&[("a", "ACGT")]));
}

#[test]
fn checks_the_whole_input() {
    // as vsearch: invalid input after the first records still fails, and
    // the records before are written
    let input = ">a\nACGT\n>b\nACGT\n>c\nAC-T\n";
    let output = run(&["1"], input);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a", "ACGT")])
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character '-' on line 6 of FASTA file"),
        "{stderr}"
    );
}

#[test]
fn warns_about_stripped_characters_after_the_first_records() {
    let output = run(&["1"], ">a\nACGT\n>b\nAC1T\n");
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a", "ACGT")])
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn rejects_invalid_counts() {
    // vsearch: the argument to --topn must be greater than zero
    for args in [&[][..], &["0"], &["-1"], &["x"], &["1", "2"]] {
        assert!(
            !run(args, &records()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn accepts_empty_input() {
    assert_eq!(topn(&["3"], ""), "");
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, records()).unwrap();
    let stdout = topn(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "1",
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(written, fasta(&[("a;size=3 desc", "ACGT")]));
}

/// Runs a chain of rsearch binaries on `input`, each with its arguments,
/// and returns the last one's output.
fn chain(steps: &[(&str, &[&str])], input: &str) -> String {
    steps.iter().fold(input.to_string(), |text, (bin, args)| {
        stdout_of(bin, run_bin(bin, args, &text))
    })
}

#[test]
fn pipeline_dereplication() {
    // vsearch --fastx_uniques in --sizein --sizeout --notrunclabels
    // --minseqlength 0 --topn 2
    let input = fasta(&[
        ("a", "ACGT"),
        ("b", "ACGT"),
        ("c", "TT"),
        ("d", "GGG"),
        ("e", "GGG"),
        ("f;size=3", "CCCCC"),
    ]);
    let derep = env!("CARGO_BIN_EXE_derep");
    let topn = env!("CARGO_BIN_EXE_topn");
    let expected = fasta(&[("f;size=3", "CCCCC"), ("a;size=2", "ACGT")]);
    assert_eq!(chain(&[(derep, &[]), (topn, &["2"])], &input), expected);
    // with --topn 100, every sequence
    let expected = fasta(&[
        ("f;size=3", "CCCCC"),
        ("a;size=2", "ACGT"),
        ("d;size=2", "GGG"),
        ("c;size=1", "TT"),
    ]);
    assert_eq!(chain(&[(derep, &[]), (topn, &["100"])], &input), expected);
}

#[test]
fn pipeline_dereplication_of_fastq() {
    // vsearch --fastx_uniques in --sizein --sizeout --notrunclabels
    // --minseqlength 0 --fastqout - --topn 2
    let input = fastq(&[
        ("a", "ACGT", "II5I"),
        ("b", "ACGT", "#III"),
        ("c", "TT", "II"),
        ("d;size=2 x", "GG", "55"),
    ]);
    let derep = env!("CARGO_BIN_EXE_derep");
    let topn = env!("CARGO_BIN_EXE_topn");
    let expected = fastq(&[("a;size=2", "ACGT", "&I7I"), ("c;size=1", "TT", "II")]);
    assert_eq!(chain(&[(derep, &[]), (topn, &["2"])], &input), expected);
}

#[test]
fn pipeline_uniques_with_filters_and_relabelling() {
    // vsearch --fastx_uniques in --minseqlength 3 --minuniquesize 2
    // --topn 1 --relabel U --sizeout
    let input = fasta(&[
        ("a", "ACGT"),
        ("b", "ACGT"),
        ("c", "TT"),
        ("d", "TT"),
        ("e", "GGG"),
        ("f", "GGG"),
        ("g", "GGG"),
    ]);
    let steps: [(&str, &[&str]); 5] = [
        (env!("CARGO_BIN_EXE_sqfilt"), &["--min-len", "3"]),
        (env!("CARGO_BIN_EXE_derep"), &["--no-sizein"]),
        (env!("CARGO_BIN_EXE_sqfilt"), &["--min-size", "2"]),
        (env!("CARGO_BIN_EXE_topn"), &["1"]),
        (
            env!("CARGO_BIN_EXE_relab"),
            &["--relabel", "U", "--sizeout"],
        ),
    ];
    assert_eq!(chain(&steps, &input), fasta(&[("U1;size=3", "GGG")]));
}
