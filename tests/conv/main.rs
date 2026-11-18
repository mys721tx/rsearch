//! End-to-end tests of the conv binary.
//!
//! Expected outputs are those of `vsearch --fastq_convert` (FASTQ input),
//! `vsearch --fasta2fastq` (FASTA input) and `vsearch --sff_convert` (SFF
//! input) on the same input, with the equivalent options.

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
fn relabels() {
    let input = fastq(&[("a;size=3 desc", "ACGT", "!+5~")]);
    let expected = fastq(&[("R1;size=1", "ACGT", "!+5~")]);
    assert_eq!(conv(&["--relabel", "R", "--sizeout"], &input), expected);
}

#[test]
fn writes_fasta_as_fastq() {
    let input = fasta(&[("a;size=3 desc", "ACGT"), ("b", "AC")]);
    let expected = fastq(&[("a;size=3 desc", "ACGT", "JJJJ"), ("b", "AC", "JJ")]);
    assert_eq!(conv(&[], &input), expected);
    // the highest output score, written with the output offset
    let expected = fastq(&[
        ("a;size=3 desc;size=1", "ACGT", "rrrr"),
        ("b;size=1", "AC", "rr"),
    ]);
    let args = ["--ascii-out", "64", "--qmax-out", "50", "--sizeout"];
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
fn relabels_sff_reads() {
    let expected = fastq(&[
        ("R1;size=1;length=6", "GTACGT", "5+!&0:"),
        ("R2;size=1;length=3", "GCC", "#$%"),
    ]);
    let input = sff_file(&sff_reads(), false, 1, b".sff");
    let args = ["--clip", "--relabel", "R", "--sizeout", "--lengthout"];
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
