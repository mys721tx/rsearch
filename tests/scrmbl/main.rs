//! End-to-end tests of the scrmbl binary.
//!
//! Expected outputs are those of `vsearch --scramble --notrunclabels
//! --fasta_width 0` with `--randseed` (`--seed`) and `--scramble_kmer`
//! (`--kmer`) on the same input. Each read gets its own SplitMix64 stream
//! from the seed and its number. With a k of 1 its bases are shuffled. With
//! a k from 2 to 9 they are rearranged keeping the count of every word of
//! up to k bases, case included. A read of k bases or fewer is left as it
//! is, and scores are never moved.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs scrmbl with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_scrmbl"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("scrmbl should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

/// Runs scrmbl and returns its standard output, which must succeed.
fn scrmbl(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "scrmbl failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-scrmbl-{}-{name}", std::process::id()))
}

const LABELS: [&str; 6] = ["s0", "s1", "s2", "s3", "low", "rep"];

/// Four random reads of 25 to 100 bases, one of mixed case and Ns, and one
/// of a repeated word.
const READS: [&str; 6] = [
    "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTG",
    "AGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCA",
    "CGTCAATACGGTTCAATGCCCTACT",
    "GCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTA",
    "acgtACGTacgtNNACGTnnacgtAC",
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
];

fn reads() -> String {
    let records: Vec<_> = LABELS.iter().copied().zip(READS).collect();
    fasta(&records)
}

/// The reads with the given sequences.
fn scrambled(seqs: [&str; 6]) -> String {
    let records: Vec<_> = LABELS.iter().copied().zip(seqs).collect();
    fasta(&records)
}

#[test]
fn shuffles_bases() {
    // a k of 1 by default
    let expected = scrambled([
        "GCCCCTCTCACGCAGAGATCCCGCCCTCCACTCAAGATGT",
        "GAAACGGAGTTGACAACGCGGTAAGATAGTGGGTGTGCACATGCGAGGTCGACTGGCGTG",
        "AGATAGGCTACTACTACTCGCCCTT",
        "GTACTGGTGGCGCTGACTCTACTGGCGCATGGATGGTTTCCTTCCTGGCCTGGAACGGTATGGAAGGCGACGAGTATTGTGCCATTGCTCGCCGCTGAGG",
        "gGTcctncaAATgaNtAnCgGCCtNa",
        "AAGGATTCCCGTGAGTGACTCCCGAAGTATACGCTATTCG",
    ]);
    assert_eq!(scrmbl(&["--seed", "1"], &reads()), expected);
    assert_eq!(scrmbl(&["--seed", "1", "--kmer", "1"], &reads()), expected);
}

#[test]
fn keeps_dinucleotide_counts() {
    for (seed, seqs) in [
        (
            "1",
            [
                "GATCCCAACCTGCACCACTCACTAGCACCTCCCTCGTGGG",
                "AGACTACCTATCCGATCAGAGTGGATTAGCAGGTGGGTAGCGGCGAGGACAGGGGTGGTA",
                "CAACATATTCGGTGCGTCCCTAACT",
                "GGGTGTTCTCACTGTATGGATTGTGCACTGCCCACTGGATCAGTCATGGCTATGGTCAGGATTGCAGTGGTCAGGGTGATCTGGCGCCCCCCGGCTGGGA",
                "acgtACGTnnacgtNNACGTacgtAC",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ],
        ),
        (
            "42",
            [
                "GCAGGCTCCCAACCTCACTGACGTCCCTATCACCACCTGG",
                "AGACTTCTGGGGTGAGATCGGAGGGGGGTATACATAGAGTGCGAGTAGCCACCAGCGGTA",
                "CCACTAACGTCCGCAATGGTATCTT",
                "GGGTGGGGCTCATGCTGTTGGTCCCTGGTCCAGTATGCACAGCATGGGACAGCCGTGGGTGTGAGATGCTCGGGTATGATCCTGGACCATTCCTGCTTCA",
                "acgtACGTnnacgtNNACGTacgtAC",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ],
        ),
        (
            "-1",
            [
                "GTACACTGCACCTCCCACTGGGAGCAATCTCCTCACCCCG",
                "AGGTCGAGGGGGGTGTTACACGGGCCAGGAGGAGCCAGTCTACGAGATATGCTGATAGTA",
                "CGGTCATCTGTTAACACCCGCTAAT",
                "GTGGGGATGTCCCTACGACATATGGCCAGCCCATGGTGACCAGGGTGGTCATGTCTCATCTTCCATGGTCTTGCTGGTTGCGGAGGCAGCTGGCTGGTGA",
                "acgtACGTnnacgtACGTacgtNNAC",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ],
        ),
    ] {
        let output = scrmbl(&["--seed", seed, "--kmer", "2"], &reads());
        assert_eq!(output, scrambled(seqs), "{seed}");
    }
}

#[test]
fn keeps_longer_word_counts() {
    for (k, seed, seqs) in [
        (
            "3",
            "42",
            [
                "GGATCACAGTCCCCCACCCGGCCTCTCAACTGCTACACTG",
                "AGGAGAGTAGGATGGGTGCGGCTTCACTGTATATACCAGGGAGAGTACGGTCCGAGGGCA",
                "CGGTTCAATACGTCAATGCCCTACT",
                "GCTTGCATCTCTGGGAGGGGGGTGTTCCCGAGATCATCCAGCACCGTGCTGGCATTGTGGCTCATGGTGACCCTCTGCCATGGTGTACATGGGAGTGGTA",
                "acgtACGTnnacgtNNACGTacgtAC",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ],
        ),
        (
            "4",
            "-1",
            [
                "GGATCACAGTCTACACTCCAACCCGGCCCCCTGCTCACTG",
                "AGTAGGGTGCTTCAGAGGAGAGTCCGAGGGTATGTATACGGAGGATACCACTGGGCGGCA",
                "CGTCAATGCCCTACGGTTCAATACT",
                "GCATCCATGGCCCACAGCTGGGTGGGGGTGTGCATCTGGACCGTGGAGTGATCTTGTACCTCATGGTGGGTTCATGCTCTGCTCCCGAGAGGGCATTGTA",
                "acgtACGTnnacgtNNACGTacgtAC",
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ],
        ),
        (
            "6",
            "42",
            [
                READS[0],
                READS[1],
                READS[2],
                "GCATGGAGAGGGTGCTGGCCCGTGATCTGCATGGGTGGGGGTGGGCATGCTCTTGTGGTTCATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTA",
                READS[4],
                READS[5],
            ],
        ),
        ("9", "1", READS),
    ] {
        let output = scrmbl(&["--seed", seed, "--kmer", k], &reads());
        assert_eq!(output, scrambled(seqs), "k {k} seed {seed}");
    }
}

#[test]
fn leaves_short_reads() {
    // with a k of 2 or more, a read of k bases or fewer has one arrangement
    let input = fasta(&[("d", "AC"), ("e", ""), ("f", "ACGT")]);
    for k in ["2", "4", "9"] {
        assert_eq!(scrmbl(&["--seed", "1", "--kmer", k], &input), input, "{k}");
    }
}

#[test]
fn keeps_scores_in_place() {
    let input = fastq(&[
        ("q;size=2", "ACGTTGCAAC", "!\"#$%&'()*"),
        ("r", "GGCCAATT", "IIII5555"),
    ]);
    let expected = fastq(&[
        ("q;size=2", "AACACGTTGC", "!\"#$%&'()*"),
        ("r", "GGCCAATT", "IIII5555"),
    ]);
    assert_eq!(scrmbl(&["--seed", "7", "--kmer", "2"], &input), expected);
    let expected = fastq(&[
        ("q;size=2", "TACCTGACAG", "!\"#$%&'()*"),
        ("r", "ACGTGTAC", "IIII5555"),
    ]);
    assert_eq!(scrmbl(&["--seed", "7"], &input), expected);
}

/// The sorted bases of each read of FASTA text.
fn compositions(text: &str) -> Vec<Vec<u8>> {
    text.lines()
        .filter(|l| !l.starts_with('>'))
        .map(|l| {
            let mut bases = l.as_bytes().to_vec();
            bases.sort_unstable();
            bases
        })
        .collect()
}

#[test]
fn scrambles_with_a_system_seed_without_one() {
    for args in [&[][..], &["--seed", "0"]] {
        let output = scrmbl(args, &reads());
        assert_eq!(compositions(&output), compositions(&reads()), "{args:?}");
    }
}

#[test]
fn rejects_invalid_arguments() {
    // vsearch: a k from 1 to 9
    for args in [
        &["--kmer", "0"][..],
        &["--kmer", "10"],
        &["--kmer", "x"],
        &["--seed", "9223372036854775808"],
    ] {
        assert!(
            !run(args, &reads()).status.success(),
            "{args:?} should fail"
        );
    }
}

#[test]
fn rejects_invalid_abundances() {
    // the reads before are written
    let input = fasta(&[("a;size=1", "ACGT"), ("z;size=0", "AC")]);
    let output = run(&["--seed", "1"], &input);
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stdout),
        fasta(&[("a;size=1", "TGAC")])
    );
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
}

#[test]
fn rejects_invalid_input() {
    let output = run(&["--seed", "1"], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("Illegal sequence character 'I'"),
        "{stderr}"
    );
}

#[test]
fn accepts_empty_input() {
    assert_eq!(scrmbl(&["--seed", "1"], ""), "");
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, reads()).unwrap();
    let stdout = scrmbl(
        &[
            "--seed",
            "1",
            "--kmer",
            "9",
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
    assert_eq!(written, reads());
}
