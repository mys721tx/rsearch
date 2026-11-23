//! End-to-end tests of the fqchar binary.
//!
//! Expected reports are those `vsearch --fastq_chars` writes to its log
//! (`--log`), without the lines the log adds before and after. A report
//! gives the number of reads, the lowest and highest score symbols and the
//! encoding they suggest. It then gives per residue its count, share and
//! longest run, and per score symbol its share and the reads that end in a
//! run of it (`--tail`, vsearch's `--fastq_tail`).

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs fqchar with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_fqchar"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("fqchar should start");
    let mut stdin = child.stdin.take().unwrap();
    let fed = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// Runs fqchar and returns its report, which must succeed.
fn fqchar(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "fqchar failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-fqchar-{}-{name}", std::process::id()))
}

/// Reads with runs of a base, Ns of two scores, lower case, U, and reads
/// ending in runs of a score.
fn reads() -> String {
    fastq(&[
        ("a", "ACGTNNacgu", "IIII#$IIII"),
        ("b", "AAAAAC", "!!!!!!"),
        ("c", "AC", "II"),
        ("d", "GGGGGNN", "5555555"),
    ])
}

/// The report on `reads()`, before its table of score symbols.
const READS_HEAD: &str = "\
Read 4 sequences.
Qmin 33, Qmax 73, Range 41
Guess: -fastq_qmin 0 -fastq_qmax 40 -fastq_ascii 33
Guess: Original Sanger format (phred+33)

Letter          N   Freq MaxRun
------ ---------- ------ ------
     A          8  32.0%      4
     C          4  16.0%      0
     G          7  28.0%      4
     N          4  16.0%      1  Q=#..5
     T          1   4.0%      0
     U          1   4.0%      0

Char  ASCII    Freq       Tails
----  -----  ------  ----------
";

/// The table of score symbols of `reads()`, with the number of reads that
/// end in a run of each.
fn reads_report(tails: [u32; 5]) -> String {
    let rows = [
        ("!", 33, "24.0"),
        ("#", 35, "4.0"),
        ("$", 36, "4.0"),
        ("5", 53, "28.0"),
        ("I", 73, "40.0"),
    ];
    let table: String = rows
        .iter()
        .zip(tails)
        .map(|((c, ascii, freq), tail)| format!(" '{c}'  {ascii:>5}  {freq:>5}%  {tail:>10}\n"))
        .collect();
    format!("{READS_HEAD}{table}")
}

#[test]
fn reports_residues_and_scores() {
    // a run's length counts its repeats, lower case is counted as upper
    // case, and N's scores range from '#' to '5'
    assert_eq!(fqchar(&[], &reads()), reads_report([1, 0, 0, 1, 1]));
}

#[test]
fn counts_tails_of_a_given_length() {
    // reads ending in a run of at least --tail symbols, 4 by default
    assert_eq!(
        fqchar(&["--tail", "1"], &reads()),
        reads_report([1, 0, 0, 1, 2])
    );
    assert_eq!(
        fqchar(&["--tail", "7"], &reads()),
        reads_report([0, 0, 0, 1, 0])
    );
}

/// The report of one read of ACGT with the given scores, for the guesses.
fn one_read(head: &str, table: &str) -> String {
    format!(
        "Read 1 sequences.\n{head}\n\
         Letter          N   Freq MaxRun\n\
         ------ ---------- ------ ------\n\
         \x20    A          1  25.0%      0\n\
         \x20    C          1  25.0%      0\n\
         \x20    G          1  25.0%      0\n\
         \x20    T          1  25.0%      0\n\n\
         Char  ASCII    Freq       Tails\n\
         ----  -----  ------  ----------\n{table}"
    )
}

#[test]
fn guesses_illumina_1_8() {
    let expected = one_read(
        "Qmin 73, Qmax 74, Range 2\n\
         Guess: -fastq_qmin 40 -fastq_qmax 41 -fastq_ascii 33\n\
         Guess: Illumina 1.8+ format (phred+33)\n",
        " 'I'     73   50.0%           0\n 'J'     74   50.0%           0\n",
    );
    assert_eq!(fqchar(&[], &fastq(&[("a", "ACGT", "IIJJ")])), expected);
}

#[test]
fn guesses_solexa() {
    let expected = one_read(
        "Qmin 59, Qmax 104, Range 46\n\
         Guess: -fastq_qmin -5 -fastq_qmax 40 -fastq_ascii 64\n\
         Guess: Solexa format (solexa+64, not supported;\n\
         \x20      convert with --fastq_convert --fastq_solexa)\n",
        " ';'     59   25.0%           0\n '@'     64   25.0%           0\n 'h'    104   50.0%           0\n",
    );
    assert_eq!(fqchar(&[], &fastq(&[("a", "ACGT", ";@hh")])), expected);
}

#[test]
fn guesses_illumina_1_3() {
    let expected = one_read(
        "Qmin 64, Qmax 104, Range 41\n\
         Guess: -fastq_qmin 0 -fastq_qmax 40 -fastq_ascii 64\n\
         Guess: Illumina 1.3+ format (phred+64)\n",
        " '@'     64   25.0%           0\n 'A'     65   25.0%           0\n 'h'    104   50.0%           0\n",
    );
    assert_eq!(fqchar(&[], &fastq(&[("a", "ACGT", "@Ahh")])), expected);
}

#[test]
fn guesses_illumina_1_5() {
    // and an N of a single score
    let expected = "\
Read 1 sequences.
Qmin 66, Qmax 104, Range 39
Guess: -fastq_qmin 2 -fastq_qmax 40 -fastq_ascii 64
Guess: Illumina 1.5+ format (phred+64)

Letter          N   Freq MaxRun
------ ---------- ------ ------
     A          1  20.0%      0
     C          1  20.0%      0
     G          1  20.0%      0
     N          1  20.0%      0  Q=B
     T          1  20.0%      0

Char  ASCII    Freq       Tails
----  -----  ------  ----------
 'B'     66   60.0%           0
 'h'    104   40.0%           0
";
    assert_eq!(fqchar(&[], &fastq(&[("a", "ACGTN", "BBhhB")])), expected);
}

#[test]
fn reports_empty_reads() {
    // no score: the lowest and highest are 0
    let expected = "\
Read 1 sequences.
Qmin 0, Qmax 0, Range 1
Guess: -fastq_qmin -33 -fastq_qmax -33 -fastq_ascii 33
Guess: Original Sanger format (phred+33)

Letter          N   Freq MaxRun
------ ---------- ------ ------

Char  ASCII    Freq       Tails
----  -----  ------  ----------
";
    assert_eq!(fqchar(&[], &fastq(&[("a", "", "")])), expected);
}

#[test]
fn reports_empty_input() {
    assert_eq!(fqchar(&[], ""), "Read 0 sequences.\n");
}

#[test]
fn refuses_fasta() {
    // vsearch: FASTQ file expected
    let output = run(&[], ">a\nAC\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn rejects_invalid_arguments() {
    // vsearch: the argument to --fastq_tail must be greater than zero
    for args in [
        &["--tail", "0"][..],
        &["--tail", "-1"],
        &["--tail", "x"],
        &["--fasta"],
    ] {
        assert!(
            !run(args, &reads()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn rejects_invalid_input() {
    let output = run(&[], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character 'I'"),
        "{stderr}"
    );
}

#[test]
fn reads_no_abundance() {
    // vsearch --fastq_chars writes no records and reads no abundance
    let input = fastq(&[("a;size=0", "ACGT", "IIJJ")]);
    assert!(run(&[], &input).status.success());
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("report.txt"));
    std::fs::write(&input, reads()).unwrap();
    let stdout = fqchar(
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
    assert_eq!(written, reads_report([1, 0, 0, 1, 1]));
}
