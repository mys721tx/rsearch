//! End-to-end tests of the revcomp binary.
//!
//! Expected outputs are those of `vsearch --fastx_revcomp` with
//! `--fastqout` (FASTQ input) or `--fastaout --fasta_width 0` (FASTA input)
//! on the same input. Each sequence is reverse-complemented, its case kept,
//! FASTQ scores reversed and headers unchanged.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs revcomp with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_revcomp"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("revcomp should start");
    let mut stdin = child.stdin.take().unwrap();
    let fed = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// Runs revcomp and returns its standard output, which must succeed.
fn revcomp(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "revcomp failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-revcomp-{}-{name}", std::process::id()))
}

#[test]
fn complements_every_nucleotide_code() {
    // case is kept, U and u become A and a, and the IUPAC codes swap in
    // pairs; S, W and N are their own complements
    let input = fasta(&[("a x", "ACGUuRYKMBVDHSWNacgtrykmbvdhswn")]);
    let expected = fasta(&[("a x", "nwsdhbvkmryacgtNWSDHBVKMRYaACGT")]);
    assert_eq!(revcomp(&[], &input), expected);
}

#[test]
fn reverses_fastq_scores() {
    let input = fastq(&[("a;size=2 desc", "ACGTNacgtn", "!\"#$%&'()*")]);
    let expected = fastq(&[("a;size=2 desc", "nacgtNACGT", "*)('&%$#\"!")]);
    assert_eq!(revcomp(&[], &input), expected);
}

#[test]
fn keeps_empty_reads() {
    let input = fastq(&[("e", "", ""), ("a", "AC", "I#")]);
    let expected = fastq(&[("e", "", ""), ("a", "GT", "#I")]);
    assert_eq!(revcomp(&[], &input), expected);
    let input = fasta(&[("e", ""), ("a", "AC")]);
    assert_eq!(revcomp(&[], &input), fasta(&[("e", ""), ("a", "GT")]));
}

#[test]
fn twice_gives_the_input_but_for_uracil() {
    // vsearch --fastx_revcomp run twice: U comes back as T
    let input = fasta(&[("a x", "ACGUuRYKMBVDHSWNacgtrykmbvdhswn")]);
    let once = revcomp(&[], &input);
    let expected = fasta(&[("a x", "ACGTtRYKMBVDHSWNacgtrykmbvdhswn")]);
    assert_eq!(revcomp(&[], &once), expected);
}

#[test]
fn unwraps_fasta_and_strips_invalid_characters() {
    let output = run(&[], ">a\nAC\nGT\n>b\nAC1G\n");
    assert!(output.status.success());
    let expected = fasta(&[("a", "ACGT"), ("b", "CGT")]);
    assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn rejects_invalid_abundances() {
    // vsearch reads every abundance
    let output = run(&[], &fasta(&[("z;size=0", "AC")]));
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
    // the reads before are written
    let input = fasta(&[("a;size=1", "AC"), ("z;size=9223372036854775808", "AC")]);
    let output = run(&[], &input);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a;size=1", "GT")])
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("range"));
}

#[test]
fn accepts_empty_input() {
    assert_eq!(revcomp(&[], ""), "");
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
    // relabelling is relab's job, and FASTA output conv's
    for args in [&["--unknown"][..], &["--relabel", "R"], &["--fasta"]] {
        let output = run(args, &fasta(&[("a", "AC")]));
        assert!(!output.status.success(), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, fastq(&[("a", "AAC", "I#5")])).unwrap();
    let stdout = revcomp(
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
    assert_eq!(written, fastq(&[("a", "GTT", "5#I")]));
}
