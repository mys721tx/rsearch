//! End-to-end tests of the relab binary.
//!
//! Expected outputs are those of `vsearch --fastx_filter` with `--fastqout`
//! (FASTQ input) or `--fastaout` and `--fasta_width 0` (FASTA input) and the
//! equivalent label options, which write every record of these inputs.
//! Empty reads, which `--fastx_filter` discards, are compared with
//! `vsearch --fastq_convert`. The pipelines are compared with the vsearch
//! commands they replace.

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

/// Runs relab with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bin(env!("CARGO_BIN_EXE_relab"), args, input)
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

/// Runs relab and returns its standard output, which must succeed.
fn relab(args: &[&str], input: &str) -> String {
    stdout_of("relab", run(args, input))
}

/// Runs `bin` with `first` on `input`, then relab with `second` on its
/// output, and returns relab's output.
fn pipe(bin: &str, first: &[&str], second: &[&str], input: &str) -> String {
    let between = stdout_of(bin, run_bin(bin, first, input));
    relab(second, &between)
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
    std::env::temp_dir().join(format!("rsearch-relab-{}-{name}", std::process::id()))
}

/// Three reads: an abundance hidden by a description, annotations to
/// remove and a trailing `;`, and an annotation that starts the header.
const LABELS: [&str; 3] = ["a;size=3 desc", "b;ee=2.5;length=9;", "size=2;c"];
const SEQS: [&str; 3] = ["ACGT", "acgu", "NNAC"];
const QUALS: [&str; 3] = ["!+5~", "IIII", "####"];

/// The three reads as FASTQ.
fn reads_fq() -> String {
    let records: Vec<_> = (0..3).map(|i| (LABELS[i], SEQS[i], QUALS[i])).collect();
    fastq(&records)
}

/// The three reads as FASTA.
fn reads_fa() -> String {
    let records: Vec<_> = (0..3).map(|i| (LABELS[i], SEQS[i])).collect();
    fasta(&records)
}

/// The three reads with new labels, as FASTQ.
fn relabelled_fq(labels: [&str; 3]) -> String {
    let records: Vec<_> = (0..3).map(|i| (labels[i], SEQS[i], QUALS[i])).collect();
    fastq(&records)
}

/// The three reads with new labels, as FASTA.
fn relabelled_fa(labels: [&str; 3]) -> String {
    let records: Vec<_> = (0..3).map(|i| (labels[i], SEQS[i])).collect();
    fasta(&records)
}

/// Checks that relab gives the three reads `labels` with `args`, in FASTQ
/// and in FASTA.
fn assert_labels(args: &[&str], labels: [&str; 3]) {
    assert_eq!(relab(args, &reads_fq()), relabelled_fq(labels), "{args:?}");
    assert_eq!(relab(args, &reads_fa()), relabelled_fa(labels), "{args:?}");
}

#[test]
fn copies_records_unchanged() {
    // the whole header is kept, description included
    assert_eq!(relab(&[], &reads_fq()), reads_fq());
    assert_eq!(relab(&[], &reads_fa()), reads_fa());
}

#[test]
fn unwraps_fasta() {
    let input = ">a\nACGT\nAC\n>b\nGG\nTT\nA\n";
    assert_eq!(relab(&[], input), fasta(&[("a", "ACGTAC"), ("b", "GGTTA")]));
}

#[test]
fn accepts_empty_input() {
    assert_eq!(relab(&[], ""), "");
    assert_eq!(relab(&["--relabel", "R"], ""), "");
}

#[test]
fn relabels_with_prefix() {
    assert_labels(&["--relabel", "R"], ["R1", "R2", "R3"]);
    assert_labels(
        &["--relabel", "R", "--relabel-keep"],
        ["R1 a;size=3 desc", "R2 b;ee=2.5;length=9;", "R3 size=2;c"],
    );
}

#[test]
fn relabels_with_empty_prefix() {
    assert_labels(&["--relabel", ""], ["1", "2", "3"]);
}

#[test]
fn relabels_with_sequence() {
    // the sequence as it is, not normalized
    assert_labels(&["--relabel-self"], ["ACGT", "acgu", "NNAC"]);
    assert_labels(
        &["--relabel-self", "--relabel-keep"],
        [
            "ACGT a;size=3 desc",
            "acgu b;ee=2.5;length=9;",
            "NNAC size=2;c",
        ],
    );
}

#[test]
fn relabels_with_sha1() {
    // digests of the normalized sequence: case and U/T do not matter
    let a = "2108994e17f6cca9ff2352ada92b6511db076034";
    let c = "a85163cf89221b8ce7b2b5491343277eae75fa80";
    assert_labels(&["--relabel-sha1"], [a, a, c]);
}

#[test]
fn relabels_with_md5() {
    let a = "f1f8f4bf413b16ad135722aa4591043e";
    let c = "573936a0b1dac00aa7197f68eacd719d";
    assert_labels(&["--relabel-md5"], [a, a, c]);
    assert_labels(
        &["--relabel-md5", "--relabel-keep"],
        [
            &format!("{a} a;size=3 desc"),
            &format!("{a} b;ee=2.5;length=9;"),
            &format!("{c} size=2;c"),
        ],
    );
}

#[test]
fn keeps_labels_without_relabelling() {
    assert_labels(&["--relabel-keep"], LABELS);
}

#[test]
fn refuses_two_relabelings() {
    for args in [
        &["--relabel", "R", "--relabel-sha1"][..],
        &["--relabel-self", "--relabel-md5"],
        &["--relabel-sha1", "--relabel-md5"],
    ] {
        assert!(!run(args, &reads_fa()).status.success(), "{args:?}");
    }
}

#[test]
fn writes_abundance() {
    // an abundance followed by a description is not read, and counts as 1
    assert_labels(
        &["--sizeout"],
        [
            "a;size=3 desc;size=1",
            "b;ee=2.5;length=9;size=1",
            ";c;size=2",
        ],
    );
}

#[test]
fn removes_abundance() {
    assert_labels(&["--xsize"], ["a;size=3 desc", "b;ee=2.5;length=9;", ";c"]);
}

#[test]
fn writes_expected_errors() {
    let expected = relabelled_fq([
        "a;size=3 desc;ee=1.1100",
        "b;length=9;ee=0.0004000",
        "size=2;c;ee=2.5238",
    ]);
    assert_eq!(relab(&["--ee"], &reads_fq()), expected);
}

#[test]
fn formats_expected_errors() {
    // about four significant digits, at every magnitude
    let reads = [
        ("q93", "A", "~", "0.0000000005012"),
        ("q93x10", "AAAAAAAAAA", "~~~~~~~~~~", "0.000000005012"),
        ("q80", "A", "q", "0.00000001000"),
        ("q70", "A", "g", "0.0000001000"),
        ("q60", "A", "]", "0.000001000"),
        ("q57", "A", "Z", "0.000001995"),
        ("q50", "A", "S", "0.00001000"),
        ("q40", "A", "I", "0.0001000"),
        ("q30", "A", "?", "0.001000"),
        ("q20", "A", "5", "0.01000"),
        ("q10", "A", "+", "0.1000"),
        ("q0", "A", "!", "1.0000"),
        (
            "q0x20",
            "AAAAAAAAAAAAAAAAAAAA",
            "!!!!!!!!!!!!!!!!!!!!",
            "20.0000",
        ),
    ];
    let input: Vec<_> = reads.iter().map(|&(l, s, q, _)| (l, s, q)).collect();
    let labels: Vec<_> = reads
        .iter()
        .map(|&(l, _, _, ee)| format!("{l};ee={ee}"))
        .collect();
    let expected: Vec<_> = reads
        .iter()
        .zip(&labels)
        .map(|(&(_, s, q, _), l)| (l.as_str(), s, q))
        .collect();
    assert_eq!(relab(&["--ee"], &fastq(&input)), fastq(&expected));
}

#[test]
fn reads_phred64_for_expected_errors() {
    let input = fastq(&[("r1;ee=9", "ACGT", "@Jh~"), ("r2;size=2", "acgun", "hhhhh")]);
    let expected = fastq(&[
        ("r1;ee=1.1001", "ACGT", "@Jh~"),
        ("r2;size=2;ee=0.0005000", "acgun", "hhhhh"),
    ]);
    assert_eq!(relab(&["--ee", "--ascii", "64"], &input), expected);
    let expected = fastq(&[
        ("R1;size=1;ee=1.1001", "ACGT", "@Jh~"),
        ("R2;size=2;ee=0.0005000", "acgun", "hhhhh"),
    ]);
    let args = ["--ee", "--ascii", "64", "--relabel", "R", "--sizeout"];
    assert_eq!(relab(&args, &input), expected);
}

#[test]
fn refuses_scores_below_zero() {
    // Phred+33 scores read as Phred+64
    let input = fastq(&[("a", "ACGT", "!+5I")]);
    let output = run(&["--ee", "--ascii", "64"], &input);
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("below"), "{stderr}");
}

#[test]
fn refuses_expected_errors_with_fasta() {
    for args in [&["--ee"][..], &["--ee", "--ascii", "64"]] {
        assert!(!run(args, &reads_fa()).status.success(), "{args:?}");
    }
}

#[test]
fn removes_expected_errors() {
    assert_labels(&["--xee"], ["a;size=3 desc", "b;length=9;", "size=2;c"]);
}

#[test]
fn writes_length() {
    assert_labels(
        &["--lengthout"],
        [
            "a;size=3 desc;length=4",
            "b;ee=2.5;length=4",
            "size=2;c;length=4",
        ],
    );
}

#[test]
fn removes_length() {
    assert_labels(&["--xlength"], ["a;size=3 desc", "b;ee=2.5", "size=2;c"]);
}

#[test]
fn removes_all_annotations() {
    assert_labels(
        &["--xsize", "--xee", "--xlength"],
        ["a;size=3 desc", "b", ";c"],
    );
}

/// Headers that test where an annotation is recognized.
const EDGE: [&str; 10] = [
    "size=3",
    "a;size=2;size=3",
    "a;ee=1.2.3;b",
    "a;length=;size=;ee=.",
    "a;ee=1.5x;ee=2",
    "a;length=5",
    ";;",
    "a;SIZE=4",
    "a;size=+4",
    "a; size=4",
];

/// The edge headers as FASTA, each with sequence AC.
fn edge_fa() -> String {
    let records: Vec<_> = EDGE.iter().map(|&label| (label, "AC")).collect();
    fasta(&records)
}

/// The edge headers relabelled.
fn edge_labels(labels: [&str; 10]) -> String {
    let records: Vec<_> = labels.iter().map(|&label| (label, "AC")).collect();
    fasta(&records)
}

#[test]
fn finds_annotations_as_vsearch() {
    // only the first valid annotation of each name counts, and only
    // after a ';' or at the start
    let args = ["--xsize", "--xee", "--xlength"];
    let expected = edge_labels([
        "",
        "a;size=3",
        "a;b",
        "a;length=;size=",
        "a;ee=1.5x",
        "a",
        ";;",
        "a;SIZE=4",
        "a;size=+4",
        "a; size=4",
    ]);
    assert_eq!(relab(&args, &edge_fa()), expected);

    let expected = edge_labels([
        ";size=3",
        "a;size=3;size=2",
        "a;ee=1.2.3;b;size=1",
        "a;length=;size=;ee=.;size=1",
        "a;ee=1.5x;ee=2;size=1",
        "a;length=5;size=1",
        ";;size=1",
        "a;SIZE=4;size=1",
        "a;size=+4;size=1",
        "a; size=4;size=1",
    ]);
    assert_eq!(relab(&["--sizeout"], &edge_fa()), expected);
}

#[test]
fn adds_sample() {
    // a header ending with ';' takes no second one
    assert_labels(
        &["--sample", "S1"],
        [
            "a;size=3 desc;sample=S1",
            "b;ee=2.5;length=9;sample=S1",
            "size=2;c;sample=S1",
        ],
    );
}

#[test]
fn appends_suffix() {
    // the suffix is written as it is, after a trailing ';' too
    assert_labels(
        &["--label-suffix", ";x=1;"],
        [
            "a;size=3 desc;x=1;",
            "b;ee=2.5;length=9;;x=1;",
            "size=2;c;x=1;",
        ],
    );
    assert_labels(
        &["--label-suffix", "_s"],
        ["a;size=3 desc_s", "b;ee=2.5;length=9;_s", "size=2;c_s"],
    );
    assert_labels(
        &["--label-suffix", ";"],
        ["a;size=3 desc;", "b;ee=2.5;length=9;;", "size=2;c;"],
    );
}

#[test]
fn accepts_empty_sample_and_suffix() {
    let args = ["--label-suffix", "", "--sample", ""];
    let expected = edge_labels([
        "size=3;sample=",
        "a;size=2;size=3;sample=",
        "a;ee=1.2.3;b;sample=",
        "a;length=;size=;ee=.;sample=",
        "a;ee=1.5x;ee=2;sample=",
        "a;length=5;sample=",
        ";;sample=",
        "a;SIZE=4;sample=",
        "a;size=+4;sample=",
        "a; size=4;sample=",
    ]);
    assert_eq!(relab(&args, &edge_fa()), expected);
}

#[test]
fn combines_annotations_in_vsearch_order() {
    // suffix, sample, size, expected errors, length; a suffix ending with
    // ';' takes the separator of the next annotation
    let expected = relabelled_fq([
        "R1;sample=S;size=1;ee=1.1100;length=4",
        "R2;sample=S;size=1;ee=0.0004000;length=4",
        "R3;sample=S;size=2;ee=2.5238;length=4",
    ]);
    let args = [
        "--relabel",
        "R",
        "--sizeout",
        "--ee",
        "--lengthout",
        "--sample",
        "S",
        "--label-suffix",
        ";",
    ];
    assert_eq!(relab(&args, &reads_fq()), expected);

    let expected = relabelled_fq([
        "R1_s;sample=S;size=1;ee=1.1100;length=4",
        "R2_s;sample=S;size=1;ee=0.0004000;length=4",
        "R3_s;sample=S;size=2;ee=2.5238;length=4",
    ]);
    let args = [
        "--relabel",
        "R",
        "--sizeout",
        "--ee",
        "--lengthout",
        "--sample",
        "S",
        "--label-suffix",
        "_s",
    ];
    assert_eq!(relab(&args, &reads_fq()), expected);

    let expected = relabelled_fq([
        "a;size=3 desc;size=1;ee=1.1100;length=4",
        "b;size=1;ee=0.0004000;length=4",
        ";c;size=2;ee=2.5238;length=4",
    ]);
    let args = ["--sizeout", "--ee", "--lengthout"];
    assert_eq!(relab(&args, &reads_fq()), expected);
}

#[test]
fn combines_relabelling_and_annotations() {
    let (a, c) = (
        "f1f8f4bf413b16ad135722aa4591043e",
        "573936a0b1dac00aa7197f68eacd719d",
    );
    // the kept label is the input's, annotations included
    assert_labels(
        &["--sizeout", "--relabel-md5", "--relabel-keep"],
        [
            &format!("{a};size=1 a;size=3 desc"),
            &format!("{a};size=1 b;ee=2.5;length=9;"),
            &format!("{c};size=2 size=2;c"),
        ],
    );
    let (a, c) = (
        "2108994e17f6cca9ff2352ada92b6511db076034",
        "a85163cf89221b8ce7b2b5491343277eae75fa80",
    );
    assert_labels(
        &["--relabel-sha1", "--sizeout", "--lengthout"],
        [
            &format!("{a};size=1;length=4"),
            &format!("{a};size=1;length=4"),
            &format!("{c};size=2;length=4"),
        ],
    );
    assert_labels(
        &["--relabel-self", "--sizeout"],
        ["ACGT;size=1", "acgu;size=1", "NNAC;size=2"],
    );
    assert_labels(
        &["--sample", "S1", "--sizeout", "--xlength"],
        [
            "a;size=3 desc;sample=S1;size=1",
            "b;ee=2.5;sample=S1;size=1",
            ";c;sample=S1;size=2",
        ],
    );
    // removing is ignored when relabelling
    assert_labels(
        &["--relabel", "R", "--xsize", "--relabel-keep"],
        ["R1 a;size=3 desc", "R2 b;ee=2.5;length=9;", "R3 size=2;c"],
    );
}

#[test]
fn relabels_empty_reads() {
    // as vsearch --fastq_convert, which keeps them
    let input = fastq(&[("e1", "", ""), ("e2;size=4", "AC", "II")]);
    let relabelled = |labels: [&str; 2]| fastq(&[(labels[0], "", ""), (labels[1], "AC", "II")]);
    assert_eq!(relab(&["--relabel", "R"], &input), relabelled(["R1", "R2"]));
    assert_eq!(
        relab(&["--relabel-sha1", "--lengthout"], &input),
        relabelled([
            "da39a3ee5e6b4b0d3255bfef95601890afd80709;length=0",
            "b1fb3bec6fdb22e19a94fe4c6c4481ccba2ee9f0;length=2",
        ])
    );
    assert_eq!(
        relab(&["--relabel-self", "--relabel-keep"], &input),
        relabelled([" e1", "AC e2;size=4"])
    );
    assert_eq!(
        relab(&["--relabel-md5", "--sizeout"], &input),
        relabelled([
            "d41d8cd98f00b204e9800998ecf8427e;size=1",
            "4144e097d2fa7a491cec2a7a4322f2bc;size=4",
        ])
    );
}

#[test]
fn rejects_zero_abundance() {
    // vsearch refuses it whether or not abundances are written
    for args in [&[][..], &["--sizeout"], &["--relabel", "R"]] {
        for input in [
            fasta(&[("a;size=0", "ACGT")]),
            fastq(&[("a;size=0", "ACGT", "IIII")]),
        ] {
            let output = run(args, &input);
            assert_eq!(output.status.code(), Some(1), "{args:?} {input}");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(stderr.contains("zero"), "{stderr}");
        }
    }
}

#[test]
fn rejects_abundance_out_of_range() {
    // the largest abundance is that of a signed 64-bit integer; the reads
    // before an error are written
    let input = fasta(&[
        ("b;size=9223372036854775807", "AC"),
        ("c;size=9223372036854775808", "A"),
    ]);
    for args in [&[][..], &["--sizeout"]] {
        let output = run(args, &input);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            fasta(&[("b;size=9223372036854775807", "AC")])
        );
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("range"), "{stderr}");
    }
    let input = fasta(&[("a;size=99999999999999999999", "ACGT")]);
    let output = run(&["--sizeout"], &input);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn rejects_invalid_input() {
    // a FASTQ record without its '+' line
    let output = run(&["--relabel", "R"], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character 'I'"),
        "{stderr}"
    );
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, reads_fa()).unwrap();
    let stdout = relab(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--relabel",
            "R",
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(written, relabelled_fa(["R1", "R2", "R3"]));
}

/// Phred+64 reads for the conversion pipeline.
fn phred64_fq() -> String {
    fastq(&[("r1;size=3 x", "ACGT", "@Jh~"), ("r2", "acgun", "hhhhh")])
}

/// The Phred+64 reads converted to Phred+33, with new labels.
fn converted(labels: [&str; 2]) -> String {
    fastq(&[(labels[0], "ACGT", "!+I_"), (labels[1], "acgun", "IIIII")])
}

#[test]
fn pipeline_fastq_convert() {
    // vsearch --fastq_convert in --fastq_ascii 64 --relabel R
    let conv = env!("CARGO_BIN_EXE_conv");
    let cases: [(&[&str], [&str; 2]); 4] = [
        (&["--relabel", "R"], ["R1", "R2"]),
        (&["--relabel", "R", "--sizeout"], ["R1;size=1", "R2;size=1"]),
        (
            &["--relabel-md5", "--relabel-keep"],
            [
                "f1f8f4bf413b16ad135722aa4591043e r1;size=3 x",
                "252fe4e1c9aa67ce660443056dfa3799 r2",
            ],
        ),
        (
            &["--lengthout", "--xsize", "--sample", "S"],
            ["r1;size=3 x;sample=S;length=4", "r2;sample=S;length=5"],
        ),
    ];
    for (args, labels) in cases {
        let output = pipe(conv, &["--ascii", "64"], args, &phred64_fq());
        assert_eq!(output, converted(labels), "{args:?}");
    }
}

#[test]
fn pipeline_dereplication() {
    // vsearch --fastx_uniques in --sizein --relabel U --sizeout
    let derep = env!("CARGO_BIN_EXE_derep");
    let input = fasta(&[
        ("a;size=2 x", "ACGT"),
        ("b", "acgt"),
        ("c;size=5", "AAAA"),
        ("d", "ACGT"),
        ("e;size=1", "AAAA"),
        ("f", "GG"),
    ]);
    let cases: [(&[&str], [&str; 3]); 3] = [
        (
            &["--relabel", "U", "--sizeout"],
            ["U1;size=6", "U2;size=4", "U3;size=1"],
        ),
        (
            &["--sizeout", "--lengthout"],
            [
                "c;size=6;length=4",
                "a;size=4;length=4",
                "f;size=1;length=2",
            ],
        ),
        (
            &["--relabel-md5", "--sizeout", "--sample", "S"],
            [
                "098890dde069e9abad63f19a0d9e1f32;sample=S;size=6",
                "f1f8f4bf413b16ad135722aa4591043e;sample=S;size=4",
                "86d8d92aba9ecf9bbf89f69cb3e49588;sample=S;size=1",
            ],
        ),
    ];
    for (args, labels) in cases {
        let expected = fasta(&[(labels[0], "AAAA"), (labels[1], "ACGT"), (labels[2], "GG")]);
        assert_eq!(pipe(derep, &[], args, &input), expected, "{args:?}");
    }
}
