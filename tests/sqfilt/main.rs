//! End-to-end tests of the sqfilt binary.
//!
//! Expected outputs are those of `vsearch --fastx_filter --notrunclabels`
//! with the equivalent filtering options. Kept reads go to `--fastqout` (or
//! `--fastaout --fasta_width 0`) and the others to `--fastqout_discarded`,
//! with `--reverse`, `--fastqout_rev` and `--fastqout_discarded_rev` for
//! pairs. The shares of uppercase and of non-N residues are compared with
//! `vsearch --fastx_mask --qmask soft --min_unmasked_pct` (and
//! `--hardmask`) on input already masked. The pipelines are compared with
//! the vsearch commands they replace.

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

/// Runs sqfilt with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_sqfilt"), args, input)
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

/// Runs sqfilt and returns its standard output, which must succeed.
fn sqfilt(args: &[&str], input: &str) -> String {
    stdout_of("sqfilt", run(args, input))
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
    std::env::temp_dir().join(format!("rsearch-sqfilt-{}-{name}", std::process::id()))
}

/// Writes `text` to a temporary file and returns its path.
fn temp_file(name: &str, text: &str) -> PathBuf {
    let path = temp_path(name);
    std::fs::write(&path, text).unwrap();
    path
}

/// Reads and removes a temporary file.
fn take(path: &PathBuf) -> String {
    let text = std::fs::read_to_string(path).unwrap();
    std::fs::remove_file(path).unwrap();
    text
}

/// Five reads: a long one with a score of 2 and a size hidden by a
/// description, one with an N, one with two Ns (one lower case) and a score
/// of 2, a clean one, and an empty one.
const READS: [(&str, &str, &str); 5] = [
    ("a;size=3 desc", "ACGTACGTAC", "IIIII5+#II"),
    ("b;size=1", "ACGNA", "IIIII"),
    ("c;size=5", "AnNC", "#III"),
    ("d", "ACGTACGT", "IIIIIIII"),
    ("e", "", ""),
];

/// The five reads.
fn reads() -> String {
    fastq(&READS)
}

/// The reads named by their first letters, in input order.
fn select(names: &str) -> String {
    let records: Vec<_> = READS
        .iter()
        .filter(|r| names.contains(&r.0[..1]))
        .copied()
        .collect();
    fastq(&records)
}

/// Checks the kept and discarded reads of sqfilt with `args`.
fn assert_split(args: &[&str], kept: &str, discarded: &str) {
    let path = temp_path(&format!("discarded-{}", args.join("_")));
    let mut all = args.to_vec();
    all.extend(["--discarded", path.to_str().unwrap()]);
    let output = sqfilt(&all, &reads());
    assert_eq!(output, select(kept), "kept with {args:?}");
    assert_eq!(take(&path), select(discarded), "discarded with {args:?}");
}

#[test]
fn drops_only_empty_reads_by_default() {
    // vsearch's default --fastq_minlen 1
    assert_split(&[], "abcd", "e");
}

#[test]
fn filters_by_length() {
    assert_split(&["--min-len", "5"], "abd", "ce");
    assert_split(&["--max-len", "5"], "bc", "ade");
}

#[test]
fn filters_by_ns() {
    // N and n both count
    assert_split(&["--max-ns", "0"], "ad", "bce");
    assert_split(&["--max-ns", "1"], "abd", "ce");
}

#[test]
fn filters_by_expected_errors() {
    assert_split(&["--max-ee", "0.5"], "bd", "ace");
    assert_split(&["--max-ee", "0.001"], "bd", "ace");
    assert_split(&["--max-ee-rate", "0.01"], "bd", "ace");
    assert_split(&["--max-ee-rate", "0.0001"], "b", "acde");
    assert_split(&["--max-ee-rate", "0"], "", "abcde");
}

#[test]
fn filters_by_lowest_score() {
    assert_split(&["--min-qual", "20"], "bd", "ace");
    assert_split(&["--min-qual", "2"], "abcd", "e");
    assert_split(&["--min-qual", "3"], "bd", "ace");
}

#[test]
fn filters_by_abundance() {
    // a size followed by a description is not read, and counts as 1
    assert_split(&["--min-size", "2"], "c", "abde");
    assert_split(&["--max-size", "3"], "abd", "ce");
}

/// Soft-masked reads: half, seven eighths, none, all and seven eighths in
/// upper case.
fn soft_masked() -> String {
    fasta(&[
        ("a x", "ACGTacgt"),
        ("b", "ACGTACGt"),
        ("c", "acgtacgt"),
        ("d", "NNACGTAC"),
        ("e", "NnACGTAC"),
    ])
}

/// The soft-masked reads named by their first letters.
fn soft(names: &str) -> String {
    let records = [
        ("a x", "ACGTacgt"),
        ("b", "ACGTACGt"),
        ("c", "acgtacgt"),
        ("d", "NNACGTAC"),
        ("e", "NnACGTAC"),
    ];
    let kept: Vec<_> = records
        .into_iter()
        .filter(|r| names.contains(&r.0[..1]))
        .collect();
    fasta(&kept)
}

#[test]
fn filters_by_share_of_uppercase() {
    // vsearch --fastx_mask --qmask soft --min_unmasked_pct and
    // --max_unmasked_pct: N counts as unmasked when upper case
    for (args, kept) in [
        (&["--min-upper-pct", "50"][..], "abde"),
        (&["--min-upper-pct", "50.1"], "bde"),
        (&["--max-upper-pct", "50"], "ac"),
        (&["--max-upper-pct", "87.5"], "abce"),
        (&["--min-upper-pct", "75", "--max-upper-pct", "99"], "be"),
    ] {
        assert_eq!(sqfilt(args, &soft_masked()), soft(kept), "{args:?}");
    }
}

/// Hard-masked reads: half, seven eighths, none and three quarters of
/// residues other than N.
fn hard(names: &str) -> String {
    let records = [
        ("a x", "ACGTNNNN"),
        ("b", "ACGTACGN"),
        ("c", "NNNNNNNN"),
        ("d", "NNACGTAC"),
    ];
    let kept: Vec<_> = records
        .into_iter()
        .filter(|r| names.contains(&r.0[..1]))
        .collect();
    fasta(&kept)
}

#[test]
fn filters_by_share_of_non_n() {
    // vsearch --fastx_mask --qmask soft --hardmask with --min_unmasked_pct
    // and --max_unmasked_pct
    for (args, kept) in [
        (&["--min-non-n-pct", "50"][..], "abd"),
        (&["--min-non-n-pct", "50.1"], "bd"),
        (&["--max-non-n-pct", "50"], "ac"),
        (&["--max-non-n-pct", "87.5"], "abcd"),
        (&["--min-non-n-pct", "75", "--max-non-n-pct", "99"], "bd"),
    ] {
        assert_eq!(sqfilt(args, &hard("abcd")), hard(kept), "{args:?}");
    }
}

#[test]
fn filters_fasta() {
    let input = fasta(&[("a x", "ACGTACGTAC"), ("b", "ACNNA"), ("e", "")]);
    let args = ["--min-len", "5", "--max-ns", "1"];
    assert_eq!(sqfilt(&args, &input), fasta(&[("a x", "ACGTACGTAC")]));
}

#[test]
fn refuses_score_options_with_fasta() {
    let input = fasta(&[("a", "ACGT")]);
    for args in [
        &["--max-ee", "1"][..],
        &["--max-ee-rate", "0.1"],
        &["--min-qual", "2"],
        &["--ascii", "64"],
        &["--qmin", "0"],
        &["--qmax", "41"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn checks_scores() {
    // as vsearch: scores out of range fail
    for (input, args, message) in [
        (
            "~III",
            &["--qmax", "41"][..],
            "quality score 93 above the maximum 41",
        ),
        (
            "I!II",
            &["--ascii", "64"],
            "quality score -31 below the minimum 0",
        ),
    ] {
        let input = fastq(&[("q", "ACGT", input)]);
        let output = run(args, &input);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{args:?}: {stderr}");
    }
    // Phred+64 scores
    let input = fastq(&[("a", "ACGT", "hhhh"), ("b", "ACGT", "hhJJ")]);
    let args = ["--ascii", "64", "--max-ee", "0.1"];
    assert_eq!(sqfilt(&args, &input), fastq(&[("a", "ACGT", "hhhh")]));
}

#[test]
fn rejects_invalid_abundances() {
    // vsearch reads every abundance, whether or not it filters by them
    for args in [&[][..], &["--min-size", "2"]] {
        let output = run(args, &fasta(&[("a;size=0", "ACGT")]));
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
        let input = fasta(&[("a;size=9223372036854775808", "ACGT")]);
        assert_eq!(run(args, &input).status.code(), Some(1), "{args:?}");
    }
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--min-len", "0"][..],
        &["--max-len", "0"],
        &["--max-ee", "0"],
        &["--max-ee", "-1"],
        &["--max-ee-rate", "-0.1"],
        &["--min-qual", "-1"],
        &["--min-size", "0"],
        &["--max-size", "0"],
        &["--min-upper-pct", "-1"],
        &["--max-upper-pct", "101"],
        &["--min-non-n-pct", "x"],
        // vsearch refuses a lowest share above the highest, though not
        // inverted lengths or sizes
        &["--min-upper-pct", "100", "--max-upper-pct", "95"],
        &["--min-non-n-pct", "60", "--max-non-n-pct", "50"],
        &["--ascii", "40"],
        &["--qmin", "10", "--qmax", "5"],
        &["--rev-out", "r.fq"],
        &["--rev-discarded", "r.fq"],
    ] {
        assert!(
            !run(args, &reads()).status.success(),
            "{args:?} should fail"
        );
    }
}

/// Forward and reverse reads of three pairs: the first has a reverse read
/// with low scores, the second a forward read ending with a score of 2,
/// the third a short reverse read.
fn pairs() -> (String, String) {
    let forward = fastq(&[
        ("p1/1", "ACGTACGT", "IIIIIIII"),
        ("p2/1", "ACGTAC", "IIIII#"),
        ("p3/1", "ACGTACGT", "IIIIIIII"),
    ]);
    let reverse = fastq(&[
        ("p1/2", "TTGGCCAA", "IIII####"),
        ("p2/2", "TTGGCCAA", "IIIIIIII"),
        ("p3/2", "TTGG", "IIII"),
    ]);
    (forward, reverse)
}

/// The reads of the given pairs, forward then reverse.
fn pair_reads(pairs: &[usize]) -> (String, String) {
    let (forward, reverse) = self::pairs();
    let pick = |text: &str| {
        let lines: Vec<&str> = text.lines().collect();
        pairs
            .iter()
            .map(|&p| lines[4 * (p - 1)..4 * p].join("\n") + "\n")
            .collect::<String>()
    };
    (pick(&forward), pick(&reverse))
}

/// Runs sqfilt on the pairs with `args` and returns the four outputs: kept
/// forward, kept reverse, discarded forward and discarded reverse.
fn filter_pairs(name: &str, forward: &str, reverse: &str, args: &[&str]) -> [String; 4] {
    let (fwd, rev) = (
        temp_file(&format!("{name}-f"), forward),
        temp_file(&format!("{name}-r"), reverse),
    );
    let outs = ["rev-out", "discarded", "rev-discarded"].map(|o| temp_path(&format!("{name}-{o}")));
    let mut all = vec![
        "-i",
        fwd.to_str().unwrap(),
        "--reverse",
        rev.to_str().unwrap(),
        "--rev-out",
        outs[0].to_str().unwrap(),
        "--discarded",
        outs[1].to_str().unwrap(),
        "--rev-discarded",
        outs[2].to_str().unwrap(),
    ];
    all.extend(args);
    let kept = sqfilt(&all, "");
    std::fs::remove_file(fwd).unwrap();
    std::fs::remove_file(rev).unwrap();
    [kept, take(&outs[0]), take(&outs[1]), take(&outs[2])]
}

#[test]
fn filters_pairs_together() {
    // a pair is discarded when either read is
    let (forward, reverse) = pairs();
    let [kept, kept_rev, discarded, discarded_rev] =
        filter_pairs("maxee", &forward, &reverse, &["--max-ee", "0.1"]);
    let (kf, kr) = pair_reads(&[3]);
    let (df, dr) = pair_reads(&[1, 2]);
    assert_eq!([kept, kept_rev, discarded, discarded_rev], [kf, kr, df, dr]);

    let [kept, kept_rev, discarded, discarded_rev] =
        filter_pairs("minlen", &forward, &reverse, &["--min-len", "5"]);
    let (kf, kr) = pair_reads(&[1, 2]);
    let (df, dr) = pair_reads(&[3]);
    assert_eq!([kept, kept_rev, discarded, discarded_rev], [kf, kr, df, dr]);
}

#[test]
fn filters_fasta_pairs() {
    let [kept, kept_rev, discarded, discarded_rev] = filter_pairs(
        "fasta",
        &fasta(&[("x", "ACGT")]),
        &fasta(&[("y", "AC")]),
        &["--min-len", "3"],
    );
    assert_eq!(
        [kept, kept_rev, discarded, discarded_rev],
        [
            String::new(),
            String::new(),
            fasta(&[("x", "ACGT")]),
            fasta(&[("y", "AC")])
        ]
    );
}

#[test]
fn fails_on_unequal_pairs() {
    // the pairs before are written
    let (forward, reverse) = pairs();
    let short = |text: &str| text.lines().take(8).collect::<Vec<_>>().join("\n") + "\n";
    for (fwd, rev, message) in [
        (
            forward.clone(),
            short(&reverse),
            "more forward reads than reverse reads",
        ),
        (
            short(&forward),
            reverse.clone(),
            "more reverse reads than forward reads",
        ),
    ] {
        let (f, r) = (temp_file("unequal-f", &fwd), temp_file("unequal-r", &rev));
        let args = ["-i", f.to_str().unwrap(), "--reverse", r.to_str().unwrap()];
        let output = run(&args, "");
        std::fs::remove_file(f).unwrap();
        std::fs::remove_file(r).unwrap();
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            pair_reads(&[1, 2]).0
        );
        let stderr = String::from_utf8_lossy(&output.stderr).to_lowercase();
        assert!(stderr.contains(message), "{stderr}");
    }
}

#[test]
fn refuses_pairs_in_different_formats() {
    let (f, r) = (
        temp_file("formats-f", &fastq(&[("x", "ACGT", "IIII")])),
        temp_file("formats-r", &fasta(&[("y", "ACGT")])),
    );
    let args = ["-i", f.to_str().unwrap(), "--reverse", r.to_str().unwrap()];
    let output = run(&args, "");
    std::fs::remove_file(f).unwrap();
    std::fs::remove_file(r).unwrap();
    assert!(!output.status.success());
}

#[test]
fn accepts_empty_input() {
    assert_eq!(sqfilt(&[], ""), "");
    assert_eq!(sqfilt(&["--max-ee", "1"], ""), "");
}

#[test]
fn rejects_invalid_input() {
    let output = run(&["--min-len", "2"], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character 'I'"),
        "{stderr}"
    );
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_file("in.fq", &reads()), temp_path("out.fq"));
    let args = [
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--max-ns",
        "0",
    ];
    assert_eq!(sqfilt(&args, ""), "");
    std::fs::remove_file(&input).unwrap();
    assert_eq!(take(&output), select("ad"));
}

/// Runs a chain of rsearch binaries on `input`, each with its arguments,
/// and returns the last one's output.
fn chain(steps: &[(&str, &[&str])], input: &str) -> String {
    steps.iter().fold(input.to_string(), |text, (bin, args)| {
        stdout_of(bin, run_bin(bin, args, &text))
    })
}

#[test]
fn pipeline_fastx_filter() {
    // vsearch --fastx_filter in --fastq_trunclen 6 --fastq_maxee 0.5
    // --relabel R --fastaout out
    let input = fastq(&[
        ("a", "ACGTACGTAC", "IIIII5+#II"),
        ("b", "ACGTA", "IIIII"),
        ("c", "ACGTACGT", "IIII+III"),
        ("d", "AACCGGTT", "IIIIIIII"),
    ]);
    let steps: [(&str, &[&str]); 4] = [
        (env!("CARGO_BIN_EXE_sqtrim"), &["--len", "6"]),
        (
            env!("CARGO_BIN_EXE_sqfilt"),
            &["--min-len", "6", "--max-ee", "0.5"],
        ),
        (env!("CARGO_BIN_EXE_relab"), &["--relabel", "R"]),
        (env!("CARGO_BIN_EXE_conv"), &["--fasta"]),
    ];
    let expected = fasta(&[("R1", "ACGTAC"), ("R2", "ACGTAC"), ("R3", "AACCGG")]);
    assert_eq!(chain(&steps, &input), expected);
}

#[test]
fn pipeline_fastq_filter_pairs() {
    // vsearch --fastx_filter fwd --reverse rev --fastq_truncqual 2
    // --fastq_minlen 5 --fastqout out --fastqout_rev out_rev: each file
    // trimmed by sqtrim, then the pairs filtered by sqfilt
    let (forward, reverse) = pairs();
    let sqtrim = env!("CARGO_BIN_EXE_sqtrim");
    let trim = |text: &str| stdout_of("sqtrim", run_bin(sqtrim, &["--qual", "2"], text));
    let [kept, kept_rev, _, _] = filter_pairs(
        "pipeline",
        &trim(&forward),
        &trim(&reverse),
        &["--min-len", "5"],
    );
    assert_eq!(kept, fastq(&[("p2/1", "ACGTA", "IIIII")]));
    assert_eq!(kept_rev, fastq(&[("p2/2", "TTGGCCAA", "IIIIIIII")]));
}

#[test]
fn pipeline_fastx_uniques() {
    // vsearch --fastx_uniques in --minseqlength 4 --minuniquesize 2
    // --relabel U --sizeout
    let long = "ACGT".repeat(25);
    let input = fasta(&[
        ("a", "ACG"),
        ("b", "ACG"),
        ("c", &long),
        ("d", "ACGTT"),
        ("e", &long),
        ("f", "ACGTT"),
        ("g", "ACGTT"),
        ("h", "TTTT"),
    ]);
    let steps: [(&str, &[&str]); 4] = [
        (env!("CARGO_BIN_EXE_sqfilt"), &["--min-len", "4"]),
        (env!("CARGO_BIN_EXE_derep"), &[]),
        (env!("CARGO_BIN_EXE_sqfilt"), &["--min-size", "2"]),
        (
            env!("CARGO_BIN_EXE_relab"),
            &["--relabel", "U", "--sizeout"],
        ),
    ];
    let expected = fasta(&[("U1;size=3", "ACGTT"), ("U2;size=2", &long)]);
    assert_eq!(chain(&steps, &input), expected);
}
