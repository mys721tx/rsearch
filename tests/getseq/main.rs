//! End-to-end tests of the getseq binary.
//!
//! Expected outputs are those of `vsearch --fastx_getseq` (`--label`, with
//! `--label_substr_match` for `--substr`) and `vsearch --fastx_getseqs`
//! (`--labels`, `--label_word`, `--label_words`, `--label_field`) with
//! `--notrunclabels --fasta_width 0`. The records matched go to
//! `--fastaout` (getseq's output, or its `--discarded` with `--invert`),
//! the others to `--notmatched`. Labels match whole headers, ignoring ASCII
//! case, or as a substring. Words match case-sensitively between characters
//! that are not letters or digits, or between `;` with a field.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs getseq with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_getseq"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("getseq should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs getseq and returns its standard output, which must succeed.
fn getseq(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "getseq failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-getseq-{}-{name}", std::process::id()))
}

/// Writes `text` to a temporary file and returns its path.
fn temp_file(name: &str, text: &str) -> PathBuf {
    let path = temp_path(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// Records whose headers differ in case, descriptions, annotations,
/// delimiters, and one empty header.
const RECORDS: [(&str, &str); 9] = [
    ("seq1", "A"),
    ("SEQ1", "C"),
    ("seq1 desc", "G"),
    ("seq10;size=2", "T"),
    ("x;sample=s1;size=3", "AA"),
    ("x;sample=s11", "CC"),
    ("y;tag=seq1", "GG"),
    ("", "TT"),
    ("s1-2.ok", "AC"),
];

fn records() -> String {
    fasta(&RECORDS)
}

/// The records at the given positions (from 0), in input order.
fn pick(indices: &[usize]) -> String {
    let picked: Vec<_> = indices.iter().map(|&i| RECORDS[i]).collect();
    fasta(&picked)
}

/// The records not at the given positions.
fn rest(indices: &[usize]) -> String {
    let rest: Vec<_> = (0..RECORDS.len())
        .filter(|i| !indices.contains(i))
        .collect();
    pick(&rest)
}

/// Checks the records getseq selects with `args`, the rest with
/// `--invert`, and both outputs at once with `--discarded`.
fn assert_selects(args: &[&str], indices: &[usize]) {
    assert_eq!(getseq(args, &records()), pick(indices), "{args:?}");
    let mut inverted = args.to_vec();
    inverted.push("--invert");
    assert_eq!(
        getseq(&inverted, &records()),
        rest(indices),
        "{args:?} --invert"
    );
    let path = temp_path(&format!("discarded-{}", args.join("_").replace('/', "")));
    let mut both = args.to_vec();
    both.extend(["--discarded", path.to_str().unwrap()]);
    assert_eq!(getseq(&both, &records()), pick(indices), "{args:?}");
    let discarded = std::fs::read_to_string(&path).unwrap();
    std::fs::remove_file(&path).unwrap();
    assert_eq!(discarded, rest(indices), "{args:?} --discarded");
}

#[test]
fn selects_a_label_ignoring_case() {
    // the whole header: a description makes another label
    assert_selects(&["--label", "seq1"], &[0, 1]);
}

#[test]
fn selects_a_label_as_a_substring() {
    assert_selects(&["--label", "seq1", "--substr"], &[0, 1, 2, 3, 6]);
    // an empty substring is in every header but the empty one
    assert_selects(&["--label", "", "--substr"], &[0, 1, 2, 3, 4, 5, 6, 8]);
}

#[test]
fn selects_labels_from_a_file() {
    // empty lines are skipped, and case is ignored
    let labels = temp_file("labels.txt", "seq1\n\nX;SAMPLE=S11\nseq10;size=2\n");
    let path = labels.to_str().unwrap();
    assert_selects(&["--labels", path], &[0, 1, 3, 5]);
    assert_selects(&["--labels", path, "--substr"], &[0, 1, 2, 3, 5, 6]);
    std::fs::remove_file(labels).unwrap();
}

#[test]
fn keeps_carriage_returns_in_labels_files() {
    // as vsearch, only the line feed ends a label. These match nothing
    let labels = temp_file("crlf.txt", "seq1\r\nSEQ1\r\n");
    assert_selects(&["--labels", labels.to_str().unwrap()], &[]);
    std::fs::remove_file(labels).unwrap();
}

#[test]
fn selects_a_word() {
    // between characters that are not letters or digits, and case matters
    assert_selects(&["--word", "seq1"], &[0, 2, 6]);
    assert_selects(&["--word", "s1"], &[4, 8]);
    assert_selects(&["--word", "2"], &[3, 8]);
    // an empty word sits only in an empty header here
    assert_selects(&["--word", ""], &[7]);
}

#[test]
fn selects_a_word_in_a_field() {
    // name=word between semicolons
    assert_selects(&["--word", "s1", "--field", "sample"], &[4]);
    assert_selects(&["--word", "seq1", "--field", "tag"], &[6]);
}

#[test]
fn selects_words_from_a_file() {
    let words = temp_file("words.txt", "s1\nok\nseq10\n1-2\n");
    let path = words.to_str().unwrap();
    assert_selects(&["--words", path], &[3, 4, 8]);
    assert_selects(&["--words", path, "--field", "sample"], &[4]);
    std::fs::remove_file(words).unwrap();
}

#[test]
fn selects_fastq() {
    let input = fastq(&[("a", "AC", "I#"), ("b", "GG", "55"), ("A", "T", "+")]);
    let expected = fastq(&[("a", "AC", "I#"), ("A", "T", "+")]);
    assert_eq!(getseq(&["--label", "a"], &input), expected);
}

#[test]
fn needs_exactly_one_selection() {
    // vsearch: specify one label option
    let words = temp_file("one.txt", "a\n");
    let path = words.to_str().unwrap();
    for args in [
        &[][..],
        &["--label", "a", "--word", "a"],
        &["--labels", path, "--words", path],
        &["--invert"],
    ] {
        assert!(
            !run(args, &records()).status.success(),
            "{args:?} should fail"
        );
    }
    std::fs::remove_file(words).unwrap();
}

#[test]
fn refuses_options_that_do_not_apply() {
    // a field names where a word is, and substrings are of labels
    for args in [
        &["--label", "a", "--field", "x"][..],
        &["--word", "a", "--substr"],
        &["--relabel", "R"],
    ] {
        assert!(
            !run(args, &records()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn fails_on_a_missing_labels_file() {
    let output = run(&["--labels", "/nonexistent/labels.txt"], &records());
    assert_eq!(output.status.code(), Some(1));
}

#[test]
fn unwraps_fasta_and_strips_invalid_characters() {
    let output = run(&["--label", "a"], ">a\nAC\nG1T\n>b\nAA\n");
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
fn rejects_invalid_input() {
    // the records before are written
    let output = run(&["--label", "a"], ">a\nACGT\n>b\nAC-T\n");
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a", "ACGT")])
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character '-'"),
        "{stderr}"
    );
}

#[test]
fn accepts_empty_input() {
    assert_eq!(getseq(&["--label", "a"], ""), "");
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_file("in.fa", &records()), temp_path("out.fa"));
    let stdout = getseq(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--word",
            "seq1",
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(written, pick(&[0, 2, 6]));
}

#[test]
fn rejects_invalid_abundances() {
    // vsearch reads the abundance of every record, selected or not, and
    // writes the records before the one that fails
    for (bad, message) in [
        ("z;size=0", "zero"),
        ("z;size=9223372036854775808", "range"),
    ] {
        let input = fasta(&[("a;size=1", "AC"), ("b", "GG"), (bad, "TT"), ("c", "CC")]);
        let path = temp_path(&format!("abundance-{message}"));
        let args = ["--label", "b", "--discarded", path.to_str().unwrap()];
        let output = run(&args, &input);
        let discarded = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(output.status.code(), Some(1), "{bad}");
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            fasta(&[("b", "GG")])
        );
        assert_eq!(discarded, fasta(&[("a;size=1", "AC")]));
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{bad}: {stderr}");
        // whatever the selection
        let output = run(&["--word", "x"], &input);
        assert_eq!(output.status.code(), Some(1), "{bad}");
    }
}
