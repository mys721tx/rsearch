//! End-to-end tests of the joinpa binary.
//!
//! Expected outputs are those of `vsearch --fastq_join --fastqout` with
//! `--reverse` on the same inputs. Each forward read is followed by the
//! pad, then the reverse complement of its mate. The forward header is
//! whole, and the scores are those of the forward read, the pad and the
//! reversed mate. `--pad` and `--pad-qual` are vsearch's `--join_padgap`
//! and `--join_padgapq`. The pad scores default to Q40 for each pad symbol,
//! `I`, or `h` with `--ascii 64` (`--fastq_ascii 64`).

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-joinpa-{}-{name}", std::process::id()))
}

/// Builds FASTQ text from (label, sequence, quality) triples.
fn fastq(records: &[(&str, &str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq, qual)| format!("@{id}\n{seq}\n+\n{qual}\n"))
        .collect()
}

/// Runs joinpa on `forward`, fed on stdin, and `reverse`, from a file named
/// after `name`, with `args`.
fn run(name: &str, args: &[&str], forward: &str, reverse: &str) -> Output {
    let path = temp_path(name);
    std::fs::write(&path, reverse).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_joinpa"))
        .args(["--reverse", path.to_str().unwrap()])
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("joinpa should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(forward.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    std::fs::remove_file(path).unwrap();
    output
}

/// Runs joinpa and returns the joined reads, which must succeed.
fn join(name: &str, args: &[&str], forward: &str, reverse: &str) -> String {
    let output = run(name, args, forward, reverse);
    assert!(
        output.status.success(),
        "joinpa failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

/// Asserts that a run failed with `message`, after writing `written`.
fn assert_fails(output: &Output, message: &str, written: &str) {
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains(message), "{stderr}");
    assert_eq!(String::from_utf8_lossy(&output.stdout), written);
}

/// Two forward reads, the first with a description, and their mates.
fn forward() -> String {
    fastq(&[("a/1 x", "ACGTu", "ABCDE"), ("b;size=3", "AC", "II")])
}

fn reverse() -> String {
    fastq(&[("a/2", "GGcAr", "12345"), ("zz", "T", "#")])
}

/// The first pair joined with the default pad.
const FIRST: &str = "@a/1 x\nACGTuNNNNNNNNyTgCC\n+\nABCDEIIIIIIII54321\n";

#[test]
fn joins_mates_with_the_default_pad() {
    // mates pair by position, whatever their labels and lengths; the
    // reverse read is complemented keeping its case, and its scores reversed
    let expected = format!("{FIRST}@b;size=3\nACNNNNNNNNA\n+\nIIIIIIIIII#\n");
    assert_eq!(join("default", &[], &forward(), &reverse()), expected);
}

#[test]
fn pads_with_the_given_strings() {
    let expected = "\
@a/1 x\nACGTuyTgCC\n+\nABCDE54321\n@b;size=3\nACA\n+\nII#\n";
    assert_eq!(
        join(
            "empty",
            &["--pad", "", "--pad-qual", ""],
            &forward(),
            &reverse()
        ),
        expected
    );
    // any ASCII string is a pad
    let expected = "\
@a/1 x\nACGTuxxyTgCC\n+\nABCDE!~54321\n@b;size=3\nACxxA\n+\nII!~#\n";
    assert_eq!(
        join(
            "given",
            &["--pad", "xx", "--pad-qual", "!~"],
            &forward(),
            &reverse()
        ),
        expected
    );
}

#[test]
fn pads_scores_to_the_length_of_the_pad() {
    let expected = "\
@a/1 x\nACGTuNNNyTgCC\n+\nABCDEIII54321\n@b;size=3\nACNNNA\n+\nIIIII#\n";
    assert_eq!(
        join("pad", &["--pad", "NNN"], &forward(), &reverse()),
        expected
    );
    // Q40 is h at offset 64
    let expected = "\
@a/1 x\nACGTuNNNyTgCC\n+\nABCDEhhh54321\n@b;size=3\nACNNNA\n+\nIIhhh#\n";
    assert_eq!(
        join(
            "ascii",
            &["--pad", "NNN", "--ascii", "64"],
            &forward(),
            &reverse()
        ),
        expected
    );
}

#[test]
fn complements_every_symbol() {
    // U becomes A, keeping its case
    let reads = fastq(&[("a", "AUu", "III")]);
    assert_eq!(
        join("u", &[], &reads, &reads),
        "@a\nAUuNNNNNNNNaAT\n+\nIIIIIIIIIIIIII\n"
    );
    let forward = fastq(&[("a", "", "")]);
    let reverse = fastq(&[(
        "b",
        "ACGTRYKMBVDHSWNacgtrykmbvdhswn",
        "!\"#$%&'()*+,-./0123456789:;<=>",
    )]);
    assert_eq!(
        join("iupac", &["--pad", ""], &forward, &reverse),
        "@a\nnwsdhbvkmryacgtNWSDHBVKMRYACGT\n+\n>=<;:9876543210/.-,+*)('&%$#\"!\n"
    );
}

#[test]
fn keeps_whole_headers() {
    let forward = fastq(&[("a\tb c", "A", "I")]);
    let reverse = fastq(&[("r", "A", "I")]);
    assert_eq!(
        join("header", &[], &forward, &reverse),
        "@a\tb c\nANNNNNNNNT\n+\nIIIIIIIIII\n"
    );
}

#[test]
fn reads_scores_out_of_range() {
    // vsearch checks no score, and ! is read at offset 64
    let reads = fastq(&[("a", "AC", "!I")]);
    assert_eq!(
        join("range", &["--ascii", "64"], &reads, &reads),
        "@a\nACNNNNNNNNGT\n+\n!IhhhhhhhhI!\n"
    );
}

#[test]
fn checks_forward_abundances() {
    // vsearch reads the abundance of forward reads only
    let reads = fastq(&[("a", "A", "I"), ("b;size=0", "A", "I")]);
    let output = run("zero", &[], &reads, &reads);
    assert_fails(&output, "abundance", "@a\nANNNNNNNNT\n+\nIIIIIIIIII\n");
    let reads = fastq(&[("a", "A", "I")]);
    let zero = fastq(&[("a;size=0", "A", "I")]);
    assert_eq!(
        join("rev-zero", &[], &reads, &zero),
        "@a\nANNNNNNNNT\n+\nIIIIIIIIII\n"
    );
    // a ;size= that is not an abundance passes
    for header in ["a;size=2;", "a;size=x"] {
        let reads = fastq(&[(header, "A", "I")]);
        assert!(join("not-size", &[], &reads, &reads).starts_with(&format!("@{header}\n")));
    }
}

#[test]
fn fails_on_unequal_numbers_of_reads() {
    // the pairs before the missing mate are written
    let reverse = fastq(&[("a/2", "GGcAr", "12345")]);
    let output = run("more-forward", &[], &forward(), &reverse);
    assert_fails(&output, "More forward reads than reverse reads", FIRST);
    let forward = fastq(&[("a/1 x", "ACGTu", "ABCDE")]);
    let output = run("more-reverse", &[], &forward, &self::reverse());
    assert_fails(&output, "More reverse reads than forward reads", FIRST);
    let output = run("no-forward", &[], "", &self::reverse());
    assert_fails(&output, "More reverse reads than forward reads", "");
}

#[test]
fn joins_empty_input() {
    assert_eq!(join("nothing", &[], "", ""), "");
}

#[test]
fn rejects_fasta() {
    let fasta = ">a\nAC\n";
    for (forward, reverse) in [
        (fasta, &self::reverse()),
        (&self::forward(), fasta),
        ("", fasta),
    ] {
        let output = run("fasta", &[], forward, reverse);
        assert_fails(&output, "FASTQ file expected", "");
    }
}

#[test]
fn rejects_invalid_input() {
    // the pairs before the invalid read are written
    let reverse = fastq(&[("a", "A", "I"), ("b", "A", "II")]);
    let output = run("invalid", &[], &forward(), &reverse);
    assert_fails(
        &output,
        "equally long",
        "@a/1 x\nACGTuNNNNNNNNT\n+\nABCDEIIIIIIIII\n",
    );
    let output = run("invalid-forward", &[], "@a\nACGT\nIIII\n", &reverse);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--unknown"][..],
        &["--pad", "NN", "--pad-qual", "III"],
        // the default pad is 8 long
        &["--pad-qual", "III"],
        &["--pad", "NÄ"],
        &["--pad", "N", "--pad-qual", "Ä"],
        &["--ascii", "50"],
    ] {
        let output = run("arguments", args, &forward(), &reverse());
        assert_eq!(output.status.code(), Some(2), "{args:?} should fail");
        assert!(output.stdout.is_empty());
    }
    // --reverse is required
    let output = Command::new(env!("CARGO_BIN_EXE_joinpa"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn reads_and_writes_files() {
    let (input, out) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, forward()).unwrap();
    let stdout = join(
        "files",
        &["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()],
        "",
        &reverse(),
    );
    let written = std::fs::read_to_string(&out).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(out).unwrap();
    assert_eq!(stdout, "");
    assert!(written.starts_with(FIRST), "{written}");
}
