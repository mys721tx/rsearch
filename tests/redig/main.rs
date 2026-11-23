//! End-to-end tests of the redig binary.
//!
//! Expected outputs are those of `vsearch --cut --fasta_width 0` on the
//! same input. `--pattern` is `--cut_pattern`, a pattern of IUPAC codes
//! with `^` where the forward strand is cut and `_` where the reverse
//! strand is. Fragments of the forward strand go to `-o` (`--fastaout`),
//! fragments of the reverse strand to `--rev-out` (`--fastaout_rev`). The
//! sequences never cut, and their reverse complements, go to `--uncut` and
//! `--rev-uncut` (`--fastaout_discarded`, `--fastaout_discarded_rev`).

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-redig-{}-{name}", std::process::id()))
}

/// Builds FASTA text from (label, sequence) pairs.
fn fasta(records: &[(&str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq)| format!(">{id}\n{seq}\n"))
        .collect()
}

/// What a run wrote: its status and error, the fragments of each strand,
/// and the sequences never cut and their reverse complements.
struct Cut {
    output: Output,
    forward: String,
    reverse: String,
    uncut: String,
    rev_uncut: String,
}

/// Runs redig with `pattern` and `args` on `input`, fed on stdin, writing
/// every output; temporary files are named after `name`.
fn cut(name: &str, pattern: &str, args: &[&str], input: &str) -> Cut {
    let paths = ["rev", "uncut", "rev-uncut"].map(|p| temp_path(&format!("{name}-{p}")));
    let path = |i: usize| paths[i].to_str().unwrap().to_string();
    let mut all = vec![
        "--pattern".to_string(),
        pattern.to_string(),
        "--rev-out".to_string(),
        path(0),
        "--uncut".to_string(),
        path(1),
        "--rev-uncut".to_string(),
        path(2),
    ];
    all.extend(args.iter().map(|a| a.to_string()));
    let mut child = Command::new(env!("CARGO_BIN_EXE_redig"))
        .args(&all)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("redig should start");
    let mut stdin = child.stdin.take().unwrap();
    let fed = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    // an output not written reads as empty
    let read = |i: usize| std::fs::read_to_string(&paths[i]).unwrap_or_default();
    let result = Cut {
        forward: String::from_utf8(output.stdout.clone()).unwrap(),
        reverse: read(0),
        uncut: read(1),
        rev_uncut: read(2),
        output,
    };
    for path in &paths {
        let _ = std::fs::remove_file(path);
    }
    result
}

/// Runs redig, which must succeed.
fn cut_ok(name: &str, pattern: &str, args: &[&str], input: &str) -> Cut {
    let result = cut(name, pattern, args, input);
    assert!(
        result.output.status.success(),
        "redig failed: {}",
        String::from_utf8_lossy(&result.output.stderr)
    );
    result
}

/// EcoRI: G^AATTC on the forward strand, CTTAA^G on the reverse.
const ECORI: &str = "G^AATT_C";

#[test]
fn cuts_both_strands() {
    // headers are kept whole; lower case and IUPAC codes match the
    // pattern's codes they share a base with
    let input = fasta(&[
        ("s1 desc;size=3", "AAGAATTCTTTGAATTCAA"),
        ("s2", "ACGTACGT"),
        ("s3", "gaattcNNgaWttc"),
    ]);
    let result = cut_ok("ecori", ECORI, &[], &input);
    assert_eq!(
        result.forward,
        fasta(&[
            ("s1 desc;size=3", "AAG"),
            ("s1 desc;size=3", "AATTCTTTG"),
            ("s1 desc;size=3", "AATTCAA"),
            ("s3", "g"),
            ("s3", "aattcNNg"),
            ("s3", "aWttc"),
        ])
    );
    // the reverse strand's fragments, from its 3' end
    assert_eq!(
        result.reverse,
        fasta(&[
            ("s1 desc;size=3", "AATTCTT"),
            ("s1 desc;size=3", "AATTCAAAG"),
            ("s1 desc;size=3", "TTG"),
            ("s3", "aattc"),
            ("s3", "aaWtcNNg"),
            ("s3", "g"),
        ])
    );
    assert_eq!(result.uncut, fasta(&[("s2", "ACGTACGT")]));
    assert_eq!(result.rev_uncut, fasta(&[("s2", "ACGTACGT")]));
}

#[test]
fn cuts_at_the_ends_of_the_pattern() {
    let input = fasta(&[("a", "ACGAATTCGT")]);
    let result = cut_ok("ends", "^GAATTC_", &[], &input);
    assert_eq!(result.forward, fasta(&[("a", "AC"), ("a", "GAATTCGT")]));
    let result = cut_ok("ends-swapped", "_GAATTC^", &[], &input);
    assert_eq!(result.forward, fasta(&[("a", "ACGAATTC"), ("a", "GT")]));
    // a lower-case pattern is the same pattern
    let result = cut_ok("lower", "g^aatt_c", &[], &input);
    assert_eq!(result.forward, fasta(&[("a", "ACG"), ("a", "AATTCGT")]));
}

#[test]
fn cuts_at_overlapping_matches() {
    // a fragment that would be empty or negative is not written
    let result = cut_ok("overlap", "A^A_A", &[], &fasta(&[("o", "AAAAAA")]));
    assert_eq!(
        result.forward,
        fasta(&[("o", "A"), ("o", "A"), ("o", "A"), ("o", "A"), ("o", "AA")])
    );
    assert_eq!(
        result.reverse,
        fasta(&[("o", "TT"), ("o", "T"), ("o", "T"), ("o", "T"), ("o", "T")])
    );
    // the reverse cut may come before the forward cut
    let result = cut_ok("before", "G_GT^T", &[], &fasta(&[("o", "CCGGTTCC")]));
    assert_eq!(result.forward, fasta(&[("o", "CCGGT"), ("o", "TCC")]));
    assert_eq!(result.reverse, fasta(&[("o", "CGG"), ("o", "GGAAC")]));
}

#[test]
fn matches_any_base_with_n() {
    let result = cut_ok("any", "N^N_", &[], &fasta(&[("a", "ACGAATTCGT")]));
    let bases: Vec<(&str, &str)> = "ACGAATTCGT"
        .split("")
        .filter(|b| !b.is_empty())
        .map(|b| ("a", b))
        .collect();
    assert_eq!(result.forward, fasta(&bases));
    // N in the sequence matches every base of the pattern
    let result = cut_ok("n", ECORI, &[], &fasta(&[("n", "NNNNNN")]));
    assert_eq!(result.forward, fasta(&[("n", "N"), ("n", "NNNNN")]));
    assert_eq!(result.reverse, fasta(&[("n", "NNNNN"), ("n", "N")]));
}

#[test]
fn keeps_sequences_never_cut() {
    // empty, or shorter than the pattern
    let input = fasta(&[("e", ""), ("s", "GAATT")]);
    let result = cut_ok("uncut", ECORI, &[], &input);
    assert_eq!(result.forward, "");
    assert_eq!(result.reverse, "");
    assert_eq!(result.uncut, input);
    assert_eq!(result.rev_uncut, fasta(&[("e", ""), ("s", "AATTC")]));
}

#[test]
fn cuts_empty_input() {
    let result = cut_ok("nothing", ECORI, &[], "");
    assert_eq!(result.forward, "");
    assert_eq!(result.uncut, "");
}

#[test]
fn checks_abundances() {
    // the fragments before an invalid abundance are written
    let input = fasta(&[("a", "ACGAATTCGT"), ("b;size=0", "AC")]);
    let result = cut("zero", ECORI, &[], &input);
    assert_eq!(result.output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&result.output.stderr);
    assert!(stderr.contains("abundance"), "{stderr}");
    assert_eq!(result.forward, fasta(&[("a", "ACG"), ("a", "AATTCGT")]));
}

#[test]
fn rejects_fastq() {
    let result = cut("fastq", ECORI, &[], "@a\nACGAATTCGT\n+\nIIIIIIIIII\n");
    assert_eq!(result.output.status.code(), Some(1));
    assert!(result.forward.is_empty());
    let stderr = String::from_utf8_lossy(&result.output.stderr);
    assert!(stderr.contains("FASTA file expected"), "{stderr}");
}

#[test]
fn rejects_invalid_input() {
    let result = cut("invalid", ECORI, &[], ">a\nAC\x01GT\n>\x02b\nAC\n");
    assert_eq!(result.output.status.code(), Some(1));
}

#[test]
fn rejects_invalid_patterns() {
    // one ^ and one _, some IUPAC codes, and nothing else
    for pattern in [
        "GAATTC",
        "G^AATTC",
        "GAATT_C",
        "G^^AATT_C",
        "G^AA_TT_C",
        "^_",
        "G^AAXT_C",
        "G^AA-T_C",
        "G^AA TT_C",
        "G^AAÄT_C",
    ] {
        let result = cut("pattern", pattern, &[], &fasta(&[("a", "ACGAATTCGT")]));
        assert_eq!(
            result.output.status.code(),
            Some(2),
            "{pattern} should fail"
        );
        assert!(result.forward.is_empty());
    }
    // --pattern is required
    let output = Command::new(env!("CARGO_BIN_EXE_redig"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn reads_and_writes_files() {
    let (input, out) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, fasta(&[("a", "ACGAATTCGT")])).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_redig"))
        .args([
            "--pattern",
            ECORI,
            "-i",
            input.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(out).unwrap();
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(written, fasta(&[("a", "ACG"), ("a", "AATTCGT")]));
}
