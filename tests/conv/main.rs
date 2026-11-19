//! End-to-end tests of the conv binary.
//!
//! Expected outputs are those of `vsearch --fastq_convert` (FASTQ input),
//! `vsearch --fasta2fastq` (FASTA input) and `vsearch --sff_convert` (SFF
//! input) on the same input, with the equivalent options. FASTA output is
//! compared with `vsearch --fastx_subsample --sample_pct 100 --minseqlength
//! 0 --notrunclabels --fastaout`, which writes every read in input order,
//! and score checks with `vsearch --fastx_filter --fastaout`. Labels are
//! relab's job. conv therefore writes every header unchanged.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs conv with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_bytes(args, input.as_bytes())
}

/// Runs conv with `args`, feeding binary `input` on stdin.
fn run_bytes(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_conv"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("conv should start");
    child.stdin.take().unwrap().write_all(input).unwrap();
    child.wait_with_output().unwrap()
}

/// Runs conv and returns its standard output, which must succeed.
fn conv(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "conv failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-conv-{}-{name}", std::process::id()))
}

#[test]
fn copies_phred33_unchanged() {
    // the whole header is kept, description included
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    assert_eq!(conv(&[], &input), input);
}

#[test]
fn writes_phred64() {
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    let expected = fastq(&[("a;size=3 desc", "ACGT", "@JT~")]);
    let args = ["--ascii-out", "64", "--qmax-out", "62"];
    assert_eq!(conv(&args, &input), expected);
}

#[test]
fn reads_phred64() {
    let input = fastq(&[("a", "ACGT", "@Jh~")]);
    let expected = fastq(&[("a", "ACGT", "!+I_")]);
    assert_eq!(conv(&["--ascii", "64"], &input), expected);
}

#[test]
fn limits_written_scores() {
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    let expected = fastq(&[("a;size=3 desc", "ACGT", "&+5I")]);
    let args = ["--qmin-out", "5", "--qmax-out", "40"];
    assert_eq!(conv(&args, &input), expected);
}

#[test]
fn converts_solexa_to_phred() {
    // Solexa -5, -4, -1, 0, 1, 2, 10, 40
    let input = fastq(&[("s", "ACGTACGT", ";<?@ACJh")]);
    let expected = fastq(&[("s", "ACGTACGT", "\"\"$$%&+I")]);
    assert_eq!(conv(&["--ascii", "64", "--solexa"], &input), expected);
}

#[test]
fn writes_fasta_as_fastq() {
    let input = fasta(&[("a;size=3 desc", "ACGT"), ("b", "AC")]);
    let expected = fastq(&[("a;size=3 desc", "ACGT", "JJJJ"), ("b", "AC", "JJ")]);
    assert_eq!(conv(&[], &input), expected);
    // the highest output score, written with the output offset
    let expected = fastq(&[("a;size=3 desc", "ACGT", "rrrr"), ("b", "AC", "rr")]);
    let args = ["--ascii-out", "64", "--qmax-out", "50"];
    assert_eq!(conv(&args, &input), expected);
}

#[test]
fn stops_at_scores_out_of_range() {
    // the reads before are written
    let input = fastq(&[("a", "AC", "II"), ("b", "ACGT", "!+5~")]);
    let output = run(&["--qmax", "41"], &input);
    assert!(!output.status.success());
    assert_eq!(output.stdout, fastq(&[("a", "AC", "II")]).as_bytes());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("quality score 93 above the maximum 41 in read b"));
}

#[test]
fn rejects_invalid_encodings() {
    let input = fastq(&[("a", "ACGT", "IIII")]);
    for args in [
        &["--solexa"][..],
        &["--ascii", "40"],
        &["--qmin", "10", "--qmax", "5"],
        &["--ascii", "64", "--qmin", "-32"],
        &["--qmax", "94"],
        &["--qmin-out", "10", "--qmax-out", "5"],
        &["--qmin-out", "-1"],
        &["--ascii-out", "64", "--qmax-out", "63"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn rejects_input_score_options_with_fasta() {
    let input = fasta(&[("a", "ACGT")]);
    for args in [&["--ascii", "64"][..], &["--qmin", "0"], &["--solexa"]] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.fq"));
    std::fs::write(&input, fastq(&[("a", "ACGT", "@Jh~")])).unwrap();
    let stdout = conv(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--ascii",
            "64",
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert_eq!(written, fastq(&[("a", "ACGT", "!+I_")]));
}

#[test]
fn accepts_empty_input() {
    assert_eq!(conv(&[], ""), "");
}

/// An SFF read: name, bases, Phred scores, and the quality and adapter
/// clip points (1-based, 0 for none).
struct SffRead {
    name: &'static str,
    bases: &'static str,
    quals: Vec<u8>,
    clip_qual: (u16, u16),
    clip_adapter: (u16, u16),
}

/// Pads `bytes` with zeros to a multiple of 8, as SFF sections are.
fn pad8(mut bytes: Vec<u8>) -> Vec<u8> {
    bytes.resize(bytes.len().div_ceil(8) * 8, 0);
    bytes
}

/// Flows per read of the SFF files built here.
const FLOWS: u16 = 4;

/// Builds an SFF file (version 1, flowgram format 1, key TCAG) of `reads`,
/// with an index after the reads if `index`. `version` and `magic` can be
/// set wrong to build invalid files.
fn sff_file(reads: &[SffRead], index: bool, version: u32, magic: &[u8; 4]) -> Vec<u8> {
    let header_length = (31 + FLOWS + 4).div_ceil(8) * 8;
    let mut body = Vec::new();
    for read in reads {
        let n = read.bases.len();
        let mut read_header = Vec::new();
        read_header.extend(
            (16 + read.name.len() as u16)
                .div_ceil(8)
                .wrapping_mul(8)
                .to_be_bytes(),
        );
        read_header.extend((read.name.len() as u16).to_be_bytes());
        read_header.extend((n as u32).to_be_bytes());
        for clip in [
            read.clip_qual.0,
            read.clip_qual.1,
            read.clip_adapter.0,
            read.clip_adapter.1,
        ] {
            read_header.extend(clip.to_be_bytes());
        }
        read_header.extend(read.name.as_bytes());
        body.extend(pad8(read_header));
        let mut data = vec![0; 2 * usize::from(FLOWS)];
        data.extend(vec![1; n]);
        data.extend(read.bases.as_bytes());
        data.extend(&read.quals);
        body.extend(pad8(data));
    }
    let (offset, length, index_data) = if index {
        let data = [&b".mft1.00"[..], &[0; 5]].concat();
        (
            u64::from(header_length) + body.len() as u64,
            data.len() as u32,
            pad8(data),
        )
    } else {
        (0, 0, Vec::new())
    };
    let mut header = Vec::new();
    header.extend(magic);
    header.extend(version.to_be_bytes());
    header.extend(offset.to_be_bytes());
    header.extend(length.to_be_bytes());
    header.extend((reads.len() as u32).to_be_bytes());
    header.extend(header_length.to_be_bytes());
    header.extend(4u16.to_be_bytes());
    header.extend(FLOWS.to_be_bytes());
    header.push(1);
    header.extend(b"TACG");
    header.extend(b"TCAG");
    [pad8(header), body, index_data].concat()
}

/// Two reads: the first clipped by quality to bases 3..=8, with scores
/// from 0 to 99; the second clipped by adapter to bases 2..=4.
fn sff_reads() -> Vec<SffRead> {
    vec![
        SffRead {
            name: "r1",
            bases: "ACGTACGTAC",
            quals: vec![40, 30, 20, 10, 0, 5, 15, 25, 35, 99],
            clip_qual: (3, 8),
            clip_adapter: (0, 0),
        },
        SffRead {
            name: "r2",
            bases: "GGCCA",
            quals: vec![1, 2, 3, 4, 5],
            clip_qual: (0, 0),
            clip_adapter: (2, 4),
        },
    ]
}

/// Runs conv on SFF input, which must succeed, and returns the output.
fn conv_sff(args: &[&str], input: &[u8]) -> String {
    let output = run_bytes(args, input);
    assert!(
        output.status.success(),
        "conv failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn converts_sff_to_fastq() {
    // bases outside the clip points are lower case, and scores above the
    // highest output score are lowered to it
    let expected = fastq(&[
        ("r1", "acGTACGTac", "I?5+!&0:D~"),
        ("r2", "gGCCa", "\"#$%&"),
    ]);
    let input = sff_file(&sff_reads(), false, 1, b".sff");
    assert_eq!(conv_sff(&[], &input), expected);
    // the index is skipped
    let input = sff_file(&sff_reads(), true, 1, b".sff");
    assert_eq!(conv_sff(&[], &input), expected);
}

#[test]
fn clips_sff_reads() {
    let expected = fastq(&[("r1", "GTACGT", "5+!&0:"), ("r2", "GCC", "#$%")]);
    let input = sff_file(&sff_reads(), true, 1, b".sff");
    assert_eq!(conv_sff(&["--clip"], &input), expected);
}

#[test]
fn encodes_sff_scores() {
    let expected = fastq(&[("r1", "acGTACGTac", "h^TJEEOYch"), ("r2", "gGCCa", "EEEEE")]);
    let input = sff_file(&sff_reads(), false, 1, b".sff");
    let args = ["--ascii-out", "64", "--qmax-out", "40", "--qmin-out", "5"];
    assert_eq!(conv_sff(&args, &input), expected);
}

#[test]
fn lowercases_reads_with_inverted_clip_points() {
    let read = SffRead {
        name: "r1",
        bases: "ACGTACGT",
        quals: vec![30; 8],
        clip_qual: (6, 3),
        clip_adapter: (0, 0),
    };
    let input = sff_file(&[read], false, 1, b".sff");
    assert_eq!(
        conv_sff(&[], &input),
        fastq(&[("r1", "acgtacgt", "????????")])
    );
    // which cannot be clipped
    assert!(!run_bytes(&["--clip"], &input).status.success());
}

#[test]
fn accepts_empty_sff() {
    let input = sff_file(&[], false, 1, b".sff");
    assert_eq!(conv_sff(&[], &input), "");
}

#[test]
fn rejects_invalid_sff() {
    for input in [
        sff_file(&sff_reads(), false, 2, b".sff"),
        sff_file(&sff_reads(), false, 1, b".sfg"),
    ] {
        assert!(!run_bytes(&[], &input).status.success());
    }
}

#[test]
fn stops_at_truncated_sff() {
    // the reads before are written
    let mut input = sff_file(&sff_reads(), false, 1, b".sff");
    input.truncate(input.len() - 12);
    let output = run_bytes(&[], &input);
    assert!(!output.status.success());
    assert_eq!(
        output.stdout,
        fastq(&[("r1", "acGTACGTac", "I?5+!&0:D~")]).as_bytes()
    );
}

#[test]
fn rejects_input_score_options_with_sff() {
    let input = sff_file(&sff_reads(), false, 1, b".sff");
    for args in [
        &["--ascii", "64"][..],
        &["--qmin", "0"],
        &["--qmax", "41"],
        &["--solexa"],
    ] {
        assert!(
            !run_bytes(args, &input).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn clips_only_sff() {
    let input = fastq(&[("a", "ACGT", "IIII")]);
    assert!(!run(&["--clip"], &input).status.success());
}

/// Runs conv on `input`, which must fail as vsearch does: exit code 1 and
/// an error message, rather than a panic.
fn assert_fails(args: &[&str], input: &[u8]) {
    let output = run_bytes(args, input);
    assert_eq!(output.status.code(), Some(1), "{input:?} should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.starts_with("Error: "), "{input:?}: {stderr}");
}

#[test]
fn rejects_unrecognised_input() {
    // too short to be SFF, and neither FASTA nor FASTQ
    for input in [&b"."[..], b".sf", b"xyz"] {
        assert_fails(&[], input);
    }
}

#[test]
fn rejects_malformed_fastq() {
    for input in [
        // fewer scores than bases
        &b"@a\nACGT\n+\nII\n"[..],
        // no '+' line
        b"@a\nACGT\n",
    ] {
        assert_fails(&[], input);
    }
}

#[test]
fn writes_reads_before_malformed_fastq() {
    let output = run(&[], "@a\nAC\n+\nII\n@b\nACGT\n+\nII\n");
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(output.stdout, fastq(&[("a", "AC", "II")]).as_bytes());
}

/// Runs conv on `input`, which must fail with exit code 1 and an error
/// message containing `message`.
fn assert_fails_with(args: &[&str], input: &str, message: &str) {
    let output = run(args, input);
    assert_eq!(output.status.code(), Some(1), "{input:?} should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains(message), "{input:?}: {stderr}");
}

#[test]
fn strips_carriage_returns() {
    // the header ends at the carriage return
    let input = "@a b\r\nAcGt\r\n+\r\nIIII\r\n";
    assert_eq!(conv(&[], input), fastq(&[("a b", "AcGt", "IIII")]));
}

#[test]
fn rejects_illegal_fastq_sequence_characters() {
    for symbol in ["1", "*", "-", ".", " ", "\t"] {
        let input = format!("@a\nAC{symbol}T\n+\nIIII\n");
        assert_fails_with(&[], &input, "Illegal sequence character");
    }
}

#[test]
fn rejects_illegal_quality_characters() {
    for symbol in [" ", "\t", "\x7f"] {
        let input = format!("@a\nACGT\n+\nII{symbol}I\n");
        assert_fails_with(&[], &input, "Illegal quality character");
    }
}

#[test]
fn checks_the_plus_line() {
    // empty or the header again
    let input = "@a x\nACGT\n+a x\nIIII\n";
    assert_eq!(conv(&[], input), fastq(&[("a x", "ACGT", "IIII")]));
    let input = "@a x\nACGT\n+a\nIIII\n";
    assert_fails_with(&[], input, "'+' line must be empty or identical to header");
}

#[test]
fn rejects_incomplete_fastq() {
    for input in ["@a", "@a\nACGT", "@a\nACGT\n+"] {
        assert_fails_with(&[], input, "Unexpected end of file");
    }
    // the last quality line needs no line end
    let input = "@a\nACGT\n+\nIIII";
    assert_eq!(conv(&[], input), fastq(&[("a", "ACGT", "IIII")]));
}

#[test]
fn reads_multiline_fastq() {
    let input = "@a\nAC\nGT\n+\nII\nII\n@b\nA\n+\n@\n";
    let expected = fastq(&[("a", "ACGT", "IIII"), ("b", "A", "@")]);
    assert_eq!(conv(&[], input), expected);
}

#[test]
fn rejects_illegal_header_characters() {
    for input in [
        "@a\x0bb\nACGT\n+\nIIII\n",
        "@a\x7fb\nACGT\n+\nIIII\n",
        ">a\x0bb\nACGT\n",
    ] {
        assert_fails_with(&[], input, "Illegal header character");
    }
}

#[test]
fn accepts_empty_headers() {
    let input = fastq(&[("", "ACGT", "IIII")]);
    assert_eq!(conv(&[], &input), input);
    assert_eq!(conv(&[], ">\nACGT\n"), fastq(&[("", "ACGT", "JJJJ")]));
}

#[test]
fn strips_invalid_fasta_characters_with_a_warning() {
    let output = run(&[], ">a\nAC1 G*T\n>b\nA2A\n");
    assert!(output.status.success());
    let expected = fastq(&[("a", "ACGT", "JJJJ"), ("b", "AA", "JJ")]);
    assert_eq!(output.stdout, expected.as_bytes());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("4 invalid characters stripped from FASTA file:  (1) *(1) 1(1) 2(1)"),
        "{stderr}"
    );
}

#[test]
fn strips_fasta_whitespace_silently() {
    let output = run(&[], ">a\nA\tC\x0bG\x0cT\r\n");
    assert!(output.status.success());
    assert_eq!(output.stdout, fastq(&[("a", "ACGT", "JJJJ")]).as_bytes());
    assert!(output.stderr.is_empty());
}

#[test]
fn rejects_illegal_fasta_sequence_characters() {
    for symbol in ["-", ".", "\x01"] {
        let input = format!(">a\nAC{symbol}T\n");
        assert_fails_with(&[], &input, "Illegal sequence character");
    }
}

#[test]
fn requires_a_line_end_after_fasta_headers() {
    assert_fails_with(&[], ">a", "header must be terminated with newline");
    // the last sequence line needs none
    assert_eq!(conv(&[], ">a\nACGT"), fastq(&[("a", "ACGT", "JJJJ")]));
}

/// Wraps a FASTA sequence at `width`, as vsearch does: no line for an
/// empty sequence.
fn wrapped(label: &str, seq: &str, width: usize) -> String {
    let lines: String = seq
        .as_bytes()
        .chunks(width)
        .map(|line| format!("{}\n", std::str::from_utf8(line).unwrap()))
        .collect();
    format!(">{label}\n{lines}")
}

/// A long read, a short one with a description, and an empty one.
fn mixed_fq() -> String {
    fastq(&[
        ("a;size=3 desc", "ACGT", "!+5~"),
        ("e", "", ""),
        ("long", &"ACGT".repeat(45), &"I".repeat(180)),
    ])
}

#[test]
fn writes_fastq_as_fasta() {
    // vsearch --fastx_subsample --sample_pct 100 --minseqlength 0
    // --notrunclabels --fastaout, which keeps every read: lines of 80 by
    // default
    let long = "ACGT".repeat(45);
    let expected = [
        wrapped("a;size=3 desc", "ACGT", 80),
        wrapped("e", "", 80),
        wrapped("long", &long, 80),
    ]
    .concat();
    assert_eq!(conv(&["--fasta"], &mixed_fq()), expected);
    let expected = [
        wrapped("a;size=3 desc", "ACGT", 7),
        wrapped("e", "", 7),
        wrapped("long", &long, 7),
    ]
    .concat();
    assert_eq!(
        conv(&["--fasta", "--fasta-width", "7"], &mixed_fq()),
        expected
    );
}

#[test]
fn writes_unwrapped_fasta_with_width_zero() {
    // vsearch --fasta_width 0: an empty sequence has an empty line
    let expected = fasta(&[
        ("a;size=3 desc", "ACGT"),
        ("e", ""),
        ("long", &"ACGT".repeat(45)),
    ]);
    let args = ["--fasta", "--fasta-width", "0"];
    assert_eq!(conv(&args, &mixed_fq()), expected);
}

#[test]
fn rewraps_fasta() {
    let output = run(
        &["--fasta", "--fasta-width", "3"],
        ">a x\nAC\nGT\n>e\n>w\nAC1GT\n",
    );
    assert!(output.status.success());
    let expected = [
        wrapped("a x", "ACGT", 3),
        wrapped("e", "", 3),
        wrapped("w", "ACGT", 3),
    ];
    assert_eq!(String::from_utf8_lossy(&output.stdout), expected.concat());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"),
        "{stderr}"
    );
}

#[test]
fn checks_scores_when_writing_fasta() {
    // as vsearch --fastx_filter --fastaout, with --fastq_qmax 41 or
    // --fastq_ascii 64
    let input = fastq(&[("a", "ACGT", "!+5~")]);
    for (args, message) in [
        (
            &["--fasta", "--qmax", "41"][..],
            "quality score 93 above the maximum 41",
        ),
        (
            &["--fasta", "--ascii", "64"],
            "quality score -31 below the minimum 0",
        ),
    ] {
        let output = run(args, &input);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{args:?}: {stderr}");
    }
    // Phred+64 and Solexa scores are read and checked, then dropped
    let input = fastq(&[("s", "ACGT", ";@Jh")]);
    let args = ["--fasta", "--ascii", "64", "--solexa"];
    assert_eq!(conv(&args, &input), fasta(&[("s", "ACGT")]));
}

#[test]
fn refuses_output_scores_with_fasta() {
    let input = fastq(&[("a", "ACGT", "IIII")]);
    for args in [
        &["--fasta", "--ascii-out", "64"][..],
        &["--fasta", "--qmin-out", "2"],
        &["--fasta", "--qmax-out", "40"],
        &["--fasta-width", "80"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn writes_sff_as_fasta() {
    // the bases vsearch --sff_convert writes, without their scores
    let input = sff_file(&sff_reads(), true, 1, b".sff");
    let expected = fasta(&[("r1", "acGTACGTac"), ("r2", "gGCCa")]);
    assert_eq!(conv_sff(&["--fasta"], &input), expected);
    let expected = fasta(&[("r1", "GTACGT"), ("r2", "GCC")]);
    assert_eq!(conv_sff(&["--fasta", "--clip"], &input), expected);
}

#[test]
fn leaves_labels_to_relab() {
    let input = fastq(&[("a", "ACGT", "IIII")]);
    for args in [
        &["--relabel", "R"][..],
        &["--relabel-self"],
        &["--relabel-sha1"],
        &["--relabel-md5"],
        &["--relabel-keep"],
        &["--sizeout"],
        &["--xsize"],
        &["--xee"],
        &["--lengthout"],
        &["--xlength"],
        &["--sample", "S"],
        &["--label-suffix", "x"],
    ] {
        assert!(
            !run(args, &input).status.success(),
            "{args:?} should be refused"
        );
    }
}

/// Runs `bin` with `first` on `input`, then conv with `second` on its
/// output, and returns conv's output.
fn pipe(bin: &str, first: &[&str], second: &[&str], input: &str) -> String {
    let mut child = Command::new(bin)
        .args(first)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("the first binary should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let between = child.wait_with_output().unwrap();
    assert!(between.status.success(), "{bin} failed");
    conv(second, std::str::from_utf8(&between.stdout).unwrap())
}

#[test]
fn pipeline_derep_smallmem() {
    // vsearch --derep_smallmem in --sizein --sizeout --fastaout out
    let long = "ACGT".repeat(45);
    let input = fasta(&[
        ("a;size=2", &long),
        ("b", "AAAA"),
        ("c", &long),
        ("d;size=4", "AAAA"),
    ]);
    let derep = env!("CARGO_BIN_EXE_derep");
    let expected = [
        wrapped("a;size=3", &long, 80),
        wrapped("b;size=5", "AAAA", 80),
    ]
    .concat();
    assert_eq!(pipe(derep, &["--unsorted"], &["--fasta"], &input), expected);
}

#[test]
fn pipeline_fastx_filter_relabel() {
    // vsearch --fastx_filter in --relabel R --sizeout --fastaout out
    let long = "ACGT".repeat(45);
    let input = fastq(&[
        ("x;size=3 desc", "ACGT", "!+5~"),
        ("y", &long, &"I".repeat(180)),
    ]);
    let relab = env!("CARGO_BIN_EXE_relab");
    let expected = [
        wrapped("R1;size=1", "ACGT", 80),
        wrapped("R2;size=1", &long, 80),
    ]
    .concat();
    let args = ["--relabel", "R", "--sizeout"];
    assert_eq!(pipe(relab, &args, &["--fasta"], &input), expected);
}

#[test]
fn rejects_invalid_abundances() {
    // vsearch --fastq_convert and --fasta2fastq read the abundance of every
    // read, before its scores, and write the reads before the one that fails
    for (bad, message) in [
        ("z;size=0", "zero"),
        ("z;size=9223372036854775808", "range"),
    ] {
        let input = fastq(&[("a;size=1", "ACGT", "IIII"), (bad, "ACGT", "II~I")]);
        for (args, expected) in [
            (
                &["--qmax", "41"][..],
                fastq(&[("a;size=1", "ACGT", "IIII")]),
            ),
            (&["--qmax", "41", "--fasta"], fasta(&[("a;size=1", "ACGT")])),
        ] {
            let output = run(args, &input);
            assert_eq!(output.status.code(), Some(1), "{bad} {args:?}");
            assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(stderr.contains(message), "{bad} {args:?}: {stderr}");
        }
        let input = fasta(&[("a;size=1", "ACGT"), (bad, "AC")]);
        let output = run(&[], &input);
        assert_eq!(output.status.code(), Some(1), "{bad}");
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            fastq(&[("a;size=1", "ACGT", "JJJJ")])
        );
    }
}

#[test]
fn reads_no_abundance_from_sff() {
    // vsearch --sff_convert reads no abundance from read names
    let read = SffRead {
        name: "r1;size=0",
        bases: "ACGT",
        quals: vec![40, 30, 20, 10],
        clip_qual: (0, 0),
        clip_adapter: (0, 0),
    };
    let input = sff_file(&[read], false, 1, b".sff");
    assert_eq!(
        conv_sff(&[], &input),
        fastq(&[("r1;size=0", "ACGT", "I?5+")])
    );
}

// UDB input: vsearch --udb2fasta, and --fasta2fastq on its output for
// FASTQ. A UDB file holds FASTA sequences and takes FASTA input's options.
// Its headers and sequences are written as stored.

/// A sequence DUST masked, and an empty one.
const MASKED: &str = "ACGTACGTaaaaaaaaaaaaaaaaaaaaaaaaaaaaCGTACGTACGTAGCATCGATCGATCGATCG\
                      ACTAGCATCGACTAGCATCGACTGACTAGCTAGCATCGATCGACTAGCTAG";

/// The UDB file mkudb writes for `fasta`.
fn udb_file(fasta: &str) -> Vec<u8> {
    let mut child = Command::new(env!("CARGO_BIN_EXE_mkudb"))
        .args(["--word-len", "3"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("mkudb should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(fasta.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    output.stdout
}

fn masked_udb() -> Vec<u8> {
    udb_file(&fasta(&[("e", ""), ("s desc;size=4", MASKED)]))
}

/// Runs conv on a UDB file, which must succeed.
fn conv_udb(args: &[&str], input: &[u8]) -> String {
    let output = run_bytes(args, input);
    assert!(
        output.status.success(),
        "conv failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn converts_udb_to_fasta() {
    // vsearch --udb2fasta: lines of 80, no line for an empty sequence
    assert_eq!(
        conv_udb(&["--fasta"], &masked_udb()),
        format!(">e\n>s desc;size=4\n{}\n{}\n", &MASKED[..80], &MASKED[80..])
    );
    assert_eq!(
        conv_udb(&["--fasta", "--fasta-width", "0"], &masked_udb()),
        format!(">e\n\n>s desc;size=4\n{MASKED}\n")
    );
}

#[test]
fn converts_udb_to_fastq() {
    // with --fasta2fastq's score of 41
    assert_eq!(
        conv_udb(&[], &masked_udb()),
        format!(
            "@e\n\n+\n\n@s desc;size=4\n{MASKED}\n+\n{}\n",
            "J".repeat(MASKED.len())
        )
    );
    assert_eq!(
        conv_udb(&["--qmax-out", "30"], &masked_udb()),
        format!(
            "@e\n\n+\n\n@s desc;size=4\n{MASKED}\n+\n{}\n",
            "?".repeat(MASKED.len())
        )
    );
}

#[test]
fn writes_udb_headers_as_stored() {
    // vsearch reads no abundance from a UDB file, and writes control
    // characters and zeros in headers as they are
    let mut udb = masked_udb();
    let at = udb
        .windows(16)
        .position(|w| w == b"e\0s desc;size=4\0")
        .unwrap();
    udb[at] = 0x01;
    udb[at + 3] = 0;
    udb[at + 14] = b'0';
    assert_eq!(
        conv_udb(&["--fasta", "--fasta-width", "0"], &udb),
        format!(">\x01\n\n>s\0desc;size=0\n{MASKED}\n")
    );
}

#[test]
fn rejects_invalid_udb_files() {
    let udb = masked_udb();
    // vsearch reads the whole file first and writes nothing
    for (input, message) in [
        (&udb[..udb.len() - 1], "invalid UDB file: truncated"),
        (&udb[..100], "invalid UDB file: shorter than its header"),
    ] {
        let output = run_bytes(&["--fasta"], input);
        assert_eq!(output.status.code(), Some(1));
        assert!(output.stdout.is_empty());
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            format!("Error: {message}\n")
        );
    }
    let mut longer = udb.clone();
    longer.push(b'A');
    let output = run_bytes(&["--fasta"], &longer);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    // headers that are not UTF-8 are refused, as from FASTA files
    let mut udb = udb;
    let at = udb.windows(2).position(|w| w == b"e\0").unwrap();
    udb[at] = 0xe9;
    let output = run_bytes(&["--fasta"], &udb);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn rejects_input_score_options_with_udb() {
    for args in [
        &["--ascii", "64"][..],
        &["--qmin", "2"],
        &["--qmax", "40"],
        &["--solexa"],
        &["--qmin-out", "2"],
        &["--clip"],
    ] {
        let output = run_bytes(args, &masked_udb());
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn reads_udb_files() {
    let input = temp_path("in.udb");
    std::fs::write(&input, masked_udb()).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_conv"))
        .args([
            "-i",
            input.to_str().unwrap(),
            "--fasta",
            "--fasta-width",
            "0",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let _ = std::fs::remove_file(&input);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!(">e\n\n>s desc;size=4\n{MASKED}\n")
    );
}
