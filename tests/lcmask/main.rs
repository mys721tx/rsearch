//! End-to-end tests of the lcmask binary.
//!
//! Expected outputs are those of `vsearch --fastx_mask --notrunclabels
//! --minseqlength 0 --fasta_width 0` with `--fastaout` (FASTA input) or
//! `--fastqout` (FASTQ input) on the same input. Masking is DUST by default
//! (`--qmask dust`), `--hard` for `--hardmask`, and `--hard --no-dust` for
//! `--qmask soft --hardmask`. `--clear` upper-cases every residue, which is
//! what `--qmask none` means to the commands that align. Headers and FASTQ
//! scores are written unchanged.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs lcmask with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_lcmask"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("lcmask should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs lcmask and returns its standard output, which must succeed.
fn lcmask(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "lcmask failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-lcmask-{}-{name}", std::process::id()))
}

/// Short sequences: a repeat, a complex one, the shortest that can and
/// cannot be masked, an empty one, one already in mixed case with a run of
/// N, uracil, and IUPAC codes, which DUST reads as A.
fn short() -> String {
    fasta(&[
        ("a", "ACACACACACACACACACAC"),
        ("b", "ACGTTGCAAGTCCATGGACT"),
        ("c", "AAAAAAA"),
        ("d", "AAAAAAAA"),
        ("e", ""),
        ("f x", "acgtACGTacgtNNNNNNNNNNNNNNNNNNNNacgt"),
        ("g", "GATTACAGATTACAGATTACAGATTACAGATTACA"),
        ("k", "ACGUACGUACGUACGUACGU"),
        ("l", "RYKMRYKMRYKMRYKM"),
    ])
}

/// Sequences of more than one 64-residue window, each with several
/// low-complexity regions.
fn long() -> String {
    fasta(&[
        (
            "h",
            "CCGTAATGCCTTTCCCTAACAGAGTTTTTCGAACTCGTGTATATATATATATATATATATATAT\
             ATATATTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCT",
        ),
        ("i", &"A".repeat(100)),
        (
            "j",
            "TTTAGTCGTGGGATGATCAGTGGGTAAAGGCAGCAGCAGCAGCAGCAGCAGTGGCGCGGGGTAA\
             CGCGCGCTTTTTTTTTTTTTAAGGCTCAGCTGCAACGCGGAGCT",
        ),
        (
            "w",
            "TTTCCTCATGCAATTCAAAACCATGGAGAGAGAGAGAGAGAGATCCGTAATGTAGGCGAAATAT\
             TTTTTTTTTGTAAACCATTTTACGGAGGATACCAAATTCACTACTACTACTACTACTACTACTC\
             TCCTTATTCAGGACCCCCCCCCCCTAACCTGAGGT",
        ),
    ])
}

#[test]
fn soft_masks_low_complexity() {
    // the whole sequence is upper-cased first, then each region found is
    // lower-cased
    let expected = fasta(&[
        ("a", "acacacacacacacacacac"),
        ("b", "ACGTTGCAAGTCCATGGACT"),
        ("c", "AAAAAAA"),
        ("d", "aaaaaaaa"),
        ("e", ""),
        ("f x", "ACGTACGTACGTnnnnnnnnnnnnnnnnnnnnaCGT"),
        ("g", "GATTACAGATTACAGATTACAGATTACAGATTACA"),
        ("k", "ACGUACGUACGUACGUACGU"),
        ("l", "rykmrykmrykmrykm"),
    ]);
    assert_eq!(lcmask(&[], &short()), expected);
}

#[test]
fn hard_masks_low_complexity() {
    // each region found becomes N; the rest keeps its case
    let expected = fasta(&[
        ("a", "NNNNNNNNNNNNNNNNNNNN"),
        ("b", "ACGTTGCAAGTCCATGGACT"),
        ("c", "AAAAAAA"),
        ("d", "NNNNNNNN"),
        ("e", ""),
        ("f x", "acgtACGTacgtNNNNNNNNNNNNNNNNNNNNNcgt"),
        ("g", "GATTACAGATTACAGATTACAGATTACAGATTACA"),
        ("k", "ACGUACGUACGUACGUACGU"),
        ("l", "NNNNNNNNNNNNNNNN"),
    ]);
    assert_eq!(lcmask(&["--hard"], &short()), expected);
    let expected = fasta(&[("v", "acgtacgtNNNNNNNNNNNNNNNNtgca")]);
    let input = fasta(&[("v", "acgtacgtGGGGGGGGGGGGGGGGtgca")]);
    assert_eq!(lcmask(&["--hard"], &input), expected);
}

#[test]
fn masks_across_windows() {
    let expected = fasta(&[
        (
            "h",
            "CCGTAATGCCTTTCCCTAACAGAGTTTTTCGAACTCGTGtatatatatatatatatatatatat\
             atatatTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCT",
        ),
        ("i", &"a".repeat(100)),
        (
            "j",
            "TTTAGTCGTGGGATGATCAGTGGGTAAAGgcagcagcagcagcagcagcagTGGCGCGGGGTAA\
             CGCGCGCtttttttttttttAAGGCTCAGCTGCAACGCGGAGCT",
        ),
        (
            "w",
            "TTTCCTCATGCAATTCAAAACCATGgagagagagagagagagaTCCGTAATGTAGGCGAAATAt\
             tttttttttGTAAACCATTTTACGGAGGATACCAAATTCactactactactactactactactC\
             TCCTTATTCAGGAcccccccccccTAACCTGAGGT",
        ),
    ]);
    assert_eq!(lcmask(&[], &long()), expected);
    let expected = fasta(&[
        (
            "h",
            "CCGTAATGCCTTTCCCTAACAGAGTTTTTCGAACTCGTGNNNNNNNNNNNNNNNNNNNNNNNNN\
             NNNNNNTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACTGGCAGGGCT",
        ),
        ("i", &"N".repeat(100)),
        (
            "j",
            "TTTAGTCGTGGGATGATCAGTGGGTAAAGNNNNNNNNNNNNNNNNNNNNNNTGGCGCGGGGTAA\
             CGCGCGCNNNNNNNNNNNNNAAGGCTCAGCTGCAACGCGGAGCT",
        ),
        (
            "w",
            "TTTCCTCATGCAATTCAAAACCATGNNNNNNNNNNNNNNNNNNTCCGTAATGTAGGCGAAATAN\
             NNNNNNNNNGTAAACCATTTTACGGAGGATACCAAATTCNNNNNNNNNNNNNNNNNNNNNNNNC\
             TCCTTATTCAGGANNNNNNNNNNNTAACCTGAGGT",
        ),
    ]);
    assert_eq!(lcmask(&["--hard"], &long()), expected);
}

#[test]
fn hard_masks_lowercase_without_dust() {
    // only lowercase residues become N, N included
    let expected = fasta(&[
        ("a", "ACACACACACACACACACAC"),
        ("b", "ACGTTGCAAGTCCATGGACT"),
        ("c", "AAAAAAA"),
        ("d", "AAAAAAAA"),
        ("e", ""),
        ("f x", "NNNNACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNN"),
        ("g", "GATTACAGATTACAGATTACAGATTACAGATTACA"),
        ("k", "ACGUACGUACGUACGUACGU"),
        ("l", "RYKMRYKMRYKMRYKM"),
    ]);
    assert_eq!(lcmask(&["--hard", "--no-dust"], &short()), expected);
}

#[test]
fn clears_masking() {
    let input = fasta(&[("f x", "acgtACGTnnnnNNNNuU"), ("e", ""), ("r", "rykmRYKM")]);
    let expected = fasta(&[("f x", "ACGTACGTNNNNNNNNUU"), ("e", ""), ("r", "RYKMRYKM")]);
    assert_eq!(lcmask(&["--clear"], &input), expected);
}

#[test]
fn keeps_fastq_scores() {
    let input = fastq(&[
        (
            "r;size=2",
            "acacacacacacacacacacGGT",
            "IIIIIIIIIIIIIIIIIIIII#5",
        ),
        ("e", "", ""),
    ]);
    let expected = fastq(&[
        (
            "r;size=2",
            "acacacacacacacacacacGGT",
            "IIIIIIIIIIIIIIIIIIIII#5",
        ),
        ("e", "", ""),
    ]);
    assert_eq!(lcmask(&[], &input), expected);
    let expected = fastq(&[
        (
            "r;size=2",
            "NNNNNNNNNNNNNNNNNNNNGGT",
            "IIIIIIIIIIIIIIIIIIIII#5",
        ),
        ("e", "", ""),
    ]);
    assert_eq!(lcmask(&["--hard"], &input), expected);
}

#[test]
fn unwraps_fasta_and_strips_invalid_characters() {
    let output = run(&[], ">a\nACAC\nACAC\n>b\nAC1G\n");
    assert!(output.status.success());
    // eight residues of two words repeated are not enough to mask
    let expected = fasta(&[("a", "ACACACAC"), ("b", "ACG")]);
    assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn rejects_invalid_abundances() {
    // vsearch reads every abundance, whatever the masking
    for args in [&[][..], &["--hard"], &["--hard", "--no-dust"], &["--clear"]] {
        let output = run(args, &fasta(&[("z;size=0", "AC")]));
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
    }
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
fn accepts_empty_input() {
    assert_eq!(lcmask(&[], ""), "");
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
    // without DUST only hard masking does anything; clearing masks nothing;
    // filtering by the share masked is sqfilt's job
    for args in [
        &["--unknown"][..],
        &["--no-dust"],
        &["--clear", "--hard"],
        &["--clear", "--no-dust"],
        &["--min-upper-pct", "50"],
    ] {
        let output = run(args, &fasta(&[("a", "AC")]));
        assert_eq!(output.status.code(), Some(2), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, fastq(&[("a", "TTTTTTTTT", "I#5I#5I#5")])).unwrap();
    let stdout = lcmask(
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
    assert_eq!(written, fastq(&[("a", "ttttttttt", "I#5I#5I#5")]));
}
