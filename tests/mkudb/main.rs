//! End-to-end tests of the mkudb binary.
//!
//! Expected outputs are the UDB files of `vsearch --makeudb_usearch
//! --notrunclabels` on the same input, byte for byte. `--word-len` is
//! `--wordlength`. vsearch upper-cases every sequence as it reads it, then
//! masks it as `--dbmask` says. mkudb treats lower case as masked and
//! stores the sequences as given. It therefore matches `--dbmask dust` on
//! input masked by `lcmask`, and `--dbmask none` (or `soft`, the same once
//! upper-cased) on upper-case input. vsearch's length filters
//! (`--minseqlength 32`, `--maxseqlength 50000` by default) are `sqfilt`
//! before mkudb. Where a test keeps shorter sequences, vsearch was given
//! `--minseqlength 1` or `0`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

use sha1::{Digest, Sha1};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-mkudb-{}-{name}", std::process::id()))
}

/// Builds FASTA text from (label, sequence) pairs.
fn fasta(records: &[(&str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq)| format!(">{id}\n{seq}\n"))
        .collect()
}

/// Runs `bin` with `args` on `input`, fed on stdin.
fn run(bin: &str, args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(bin)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("the binary should start");
    let mut stdin = child.stdin.take().unwrap();
    let input = input.to_vec();
    // fed from a thread, so that a large input cannot block on a full pipe
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&input);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// Runs mkudb with `args` on `input`.
fn mkudb(args: &[&str], input: &str) -> Output {
    run(env!("CARGO_BIN_EXE_mkudb"), args, input.as_bytes())
}

/// Runs mkudb, which must succeed, and returns the file it wrote.
fn mkudb_ok(args: &[&str], input: &str) -> Vec<u8> {
    let output = mkudb(args, input);
    assert!(
        output.status.success(),
        "mkudb failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

/// The 2-bit code of a word: A 0, C 1, G 2, T 3.
fn code(word: &str) -> usize {
    word.bytes().fold(0, |code, base| {
        code * 4
            + match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("not a base: {base}"),
            }
    })
}

/// Builds a UDB file as vsearch writes it: a header of 50 words, the
/// number of sequences holding each of the 4^`word_len` words, the
/// sequence numbers of each word in word order, a second header, the
/// offsets of the headers, the headers each ended by a zero, the lengths of
/// the sequences, and the sequences. `words` lists each word held by a
/// sequence with the numbers of the sequences holding it.
fn udb(word_len: u32, words: &[(&str, &[u32])], records: &[(&str, &str)]) -> Vec<u8> {
    let mut out = Vec::new();
    let put = |out: &mut Vec<u8>, value: u32| out.extend(value.to_le_bytes());
    let mut header = [0u32; 50];
    header[0] = 0x5544_4246;
    header[2] = 32;
    header[4] = word_len;
    header[5] = 1;
    header[6] = 100;
    header[13] = records.len() as u32;
    header[17] = 0x0000_746e;
    header[49] = 0x5544_4266;
    for value in header {
        put(&mut out, value);
    }
    let mut lists: Vec<Vec<u32>> = vec![Vec::new(); 1 << (2 * word_len)];
    for (word, sequences) in words {
        assert_eq!(word.len(), word_len as usize);
        lists[code(word)] = sequences.to_vec();
    }
    for list in &lists {
        put(&mut out, list.len() as u32);
    }
    put(&mut out, 0x5544_4233);
    for value in lists.iter().flatten() {
        put(&mut out, *value);
    }
    let nucleotides: u64 = records.iter().map(|(_, seq)| seq.len() as u64).sum();
    let header_bytes: u64 = records.iter().map(|(h, _)| h.len() as u64 + 1).sum();
    for value in [
        0x5544_4234,
        0x005e_0db3,
        records.len() as u32,
        nucleotides as u32,
        (nucleotides >> 32) as u32,
        header_bytes as u32,
        (header_bytes >> 32) as u32,
        0x005e_0db4,
    ] {
        put(&mut out, value);
    }
    let mut offset = 0;
    for (h, _) in records {
        put(&mut out, offset);
        offset += h.len() as u32 + 1;
    }
    for (h, _) in records {
        out.extend(h.as_bytes());
        out.push(0);
    }
    for (_, seq) in records {
        put(&mut out, seq.len() as u32);
    }
    for (_, seq) in records {
        out.extend(seq.as_bytes());
    }
    out
}

/// Decodes hexadecimal text, ignoring white space.
fn unhex(text: &str) -> Vec<u8> {
    let digits: Vec<u8> = text.bytes().filter(|b| !b.is_ascii_whitespace()).collect();
    digits
        .chunks(2)
        .map(|pair| u8::from_str_radix(std::str::from_utf8(pair).unwrap(), 16).unwrap())
        .collect()
}

/// The SHA-1 digest of `bytes`, in hexadecimal.
fn sha1(bytes: &[u8]) -> String {
    Sha1::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// Upper case, U and N, an IUPAC code, a word twice in one sequence, a
/// whole header and a sequence shorter than the word length.
fn basic() -> String {
    fasta(&[
        ("s1 desc;size=2", "ACGTACGTUUGA"),
        ("s2", "NACGRTTTACGN"),
        ("s3", "AC"),
    ])
}

/// vsearch's UDB file of `basic()` with `--wordlength 3 --dbmask none
/// --minseqlength 1`.
const BASIC_UDB: &str = "
46424455000000002000000000000000030000000100000064000000000000000000000000000000000000000000000000000000030000000000000000000000
000000006e7400000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
00000000664244550000000000000000000000000000000000000000000000000200000000000000000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000
00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000
00000000010000000000000002000000000000000000000000000000000000000000000000000000010000000000000000000000000000000100000000000000
01000000020000003342445500000000010000000000000000000000000000000000000001000000000000000100000000000000000000000100000034424455
b30d5e00030000001a000000000000001500000000000000b40d5e00000000000f00000012000000733120646573633b73697a653d32007332007333000c0000
000c000000020000004143475441434754555547414e414347525454544143474e4143";

/// The words vsearch indexes for `basic()` at word length 3: U is T, and a
/// word with N or R is skipped.
const BASIC_WORDS: [(&str, &[u32]); 9] = [
    ("ACG", &[0, 1]),
    ("CGT", &[0]),
    ("GTA", &[0]),
    ("GTT", &[0]),
    ("TAC", &[0, 1]),
    ("TGA", &[0]),
    ("TTA", &[1]),
    ("TTG", &[0]),
    ("TTT", &[0, 1]),
];

#[test]
fn writes_vsearch_udb_files() {
    // the layout the other tests build their expectations with
    let records = [
        ("s1 desc;size=2", "ACGTACGTUUGA"),
        ("s2", "NACGRTTTACGN"),
        ("s3", "AC"),
    ];
    assert_eq!(udb(3, &BASIC_WORDS, &records), unhex(BASIC_UDB));
    let written = mkudb_ok(&["--word-len", "3"], &basic());
    assert_eq!(written, unhex(BASIC_UDB));
}

#[test]
fn reads_fastq() {
    let input = "@s1 desc;size=2\nACGTACGTUUGA\n+\nIIIIIIIIIIII\n\
                 @s2\nNACGRTTTACGN\n+\nIIIIIIIIIIII\n@s3\nAC\n+\nII\n";
    assert_eq!(mkudb_ok(&["--word-len", "3"], input), unhex(BASIC_UDB));
}

#[test]
fn warns_about_sequences_without_words() {
    let output = mkudb(&["--word-len", "3"], &basic());
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "Warning: 1 of 3 sequences give no word to index: they are shorter \
         than the word length (3), or each of their words holds a lower-case \
         residue or a symbol other than A, C, G, T or U\n"
    );
}

#[test]
fn skips_lowercase_words() {
    // the lcmask output of sequences vsearch masks with --dbmask dust: the
    // tandem repeat and the homopolymer; no warning, as both have words
    let records = [
        (
            "m1",
            "GCTAAAGACAATTACATAACcacacacacacacacacacacacaATACACGTCAGCAC",
        ),
        ("m2 two", "GAAACTTGTTGGttttttttttttttttCCCAGTGTGAAT"),
    ];
    let words: [(&str, &[u32]); 30] = [
        ("AAACTTGT", &[1]),
        ("AAAGACAA", &[0]),
        ("AACTTGTT", &[1]),
        ("AAGACAAT", &[0]),
        ("AATTACAT", &[0]),
        ("ACAATTAC", &[0]),
        ("ACACGTCA", &[0]),
        ("ACGTCAGC", &[0]),
        ("ACTTGTTG", &[1]),
        ("AGACAATT", &[0]),
        ("AGTGTGAA", &[1]),
        ("ATACACGT", &[0]),
        ("ATTACATA", &[0]),
        ("CAATTACA", &[0]),
        ("CACGTCAG", &[0]),
        ("CAGTGTGA", &[1]),
        ("CCAGTGTG", &[1]),
        ("CCCAGTGT", &[1]),
        ("CGTCAGCA", &[0]),
        ("CTAAAGAC", &[0]),
        ("CTTGTTGG", &[1]),
        ("GAAACTTG", &[1]),
        ("GACAATTA", &[0]),
        ("GCTAAAGA", &[0]),
        ("GTCAGCAC", &[0]),
        ("GTGTGAAT", &[1]),
        ("TAAAGACA", &[0]),
        ("TACACGTC", &[0]),
        ("TACATAAC", &[0]),
        ("TTACATAA", &[0]),
    ];
    let output = mkudb(&[], &fasta(&records));
    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert_eq!(output.stdout, udb(8, &words, &records));
}

#[test]
fn lists_sequences_in_order() {
    // twenty sequences sharing two words
    let suffixes = [
        "TTTC", "CTCA", "TGCA", "ATTC", "AAAA", "CCAT", "GTCC", "GTAA", "TGTA", "GGCG", "AAAT",
        "AGTA", "AACC", "ATTT", "TACG", "GAGG", "ATAC", "CAAA", "TTCC", "TCCT",
    ];
    let labels: Vec<String> = (0..20).map(|i| format!("r{i}")).collect();
    let seqs: Vec<String> = suffixes.iter().map(|s| format!("ACGTAC{s}")).collect();
    let records: Vec<(&str, &str)> = labels
        .iter()
        .zip(&seqs)
        .map(|(l, s)| (l.as_str(), s.as_str()))
        .collect();
    let all: Vec<u32> = (0..20).collect();
    let words: Vec<(&str, &[u32])> = vec![
        ("AAAA", &[4]),
        ("AAAT", &[10]),
        ("AACC", &[12]),
        ("ACAA", &[4, 10, 12]),
        ("ACAG", &[11]),
        ("ACAT", &[3, 13, 16]),
        ("ACCA", &[17]),
        ("ACCC", &[5]),
        ("ACCT", &[1]),
        ("ACGA", &[15]),
        ("ACGG", &[9]),
        ("ACGT", &all),
        ("ACTA", &[14]),
        ("ACTC", &[19]),
        ("ACTG", &[2, 8]),
        ("ACTT", &[0, 18]),
        ("AGTA", &[11]),
        ("ATAC", &[16]),
        ("ATTC", &[3]),
        ("ATTT", &[13]),
        ("CAAA", &[4, 10, 17]),
        ("CAAC", &[12]),
        ("CAGT", &[11]),
        ("CATA", &[16]),
        ("CATT", &[3, 13]),
        ("CCAA", &[17]),
        ("CCAT", &[5]),
        ("CCCA", &[5]),
        ("CCTC", &[1]),
        ("CGAG", &[15]),
        ("CGGC", &[9]),
        ("CGTA", &all),
        ("CGTC", &[6]),
        ("CTAC", &[14]),
        ("CTCA", &[1]),
        ("CTCC", &[19]),
        ("CTGC", &[2]),
        ("CTGT", &[8]),
        ("CTTC", &[18]),
        ("CTTT", &[0]),
        ("GAGG", &[15]),
        ("GGCG", &[9]),
        ("GTAA", &[7]),
        ("GTAC", &all),
        ("GTCC", &[6]),
        ("TACA", &[3, 4, 10, 11, 12, 13, 16]),
        ("TACC", &[1, 5, 17]),
        ("TACG", &[6, 7, 9, 14, 15]),
        ("TACT", &[0, 2, 8, 14, 18, 19]),
        ("TCCT", &[19]),
        ("TGCA", &[2]),
        ("TGTA", &[8]),
        ("TTCC", &[18]),
        ("TTTC", &[0]),
    ];
    let written = mkudb_ok(&["--word-len", "4"], &fasta(&records));
    assert_eq!(written, udb(4, &words, &records));
}

#[test]
fn keeps_empty_sequences() {
    // vsearch with --minseqlength 0; the default word length is 8
    let records = [("a", ""), ("b", "ACGTACGTAC")];
    let words: [(&str, &[u32]); 3] = [("ACGTACGT", &[1]), ("CGTACGTA", &[1]), ("GTACGTAC", &[1])];
    let output = mkudb(&[], &fasta(&records));
    assert!(output.status.success());
    assert_eq!(output.stdout, udb(8, &words, &records));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.starts_with("Warning: 1 of 2 sequences"), "{stderr}");
    assert!(stderr.contains("word length (8)"), "{stderr}");
}

#[test]
fn rejects_empty_input() {
    // vsearch writes no UDB file without sequences. None could be read back
    let output = mkudb(&[], "");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "Error: no sequences: a UDB file needs at least one\n"
    );
    // the output file is created, and left empty
    let out = temp_path("empty.udb");
    let output = mkudb(&["-o", out.to_str().unwrap()], "");
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(std::fs::read(&out).unwrap(), b"");
    let _ = std::fs::remove_file(&out);
}

#[test]
fn checks_abundances() {
    // vsearch reads every record before writing, and nothing is written
    let input = fasta(&[("a", "ACGTACGTAC"), ("b;size=0", "ACGTACGTAC")]);
    let output = mkudb(&[], &input);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("abundance"), "{stderr}");
    let input = fasta(&[("a;size=9223372036854775808", "ACGTACGTAC")]);
    assert_eq!(mkudb(&[], &input).status.code(), Some(1));
    // the largest abundance is accepted
    let input = fasta(&[("a;size=9223372036854775807", "ACGTACGTAC")]);
    assert!(mkudb(&[], &input).status.success());
}

#[test]
fn rejects_invalid_input() {
    let output = mkudb(&[], ">a\nAC\x01GT\n>\x02b\nAC\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let output = mkudb(&[], "ACGT\n");
    assert_eq!(output.status.code(), Some(1));
}

#[test]
fn rejects_invalid_word_lengths() {
    // vsearch takes 3 to 15
    for word_len in ["2", "16", "0", "-1", "x", ""] {
        let output = mkudb(&["--word-len", word_len], &basic());
        assert_eq!(output.status.code(), Some(2), "{word_len} should fail");
        assert!(output.stdout.is_empty());
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, out) = (temp_path("in.fa"), temp_path("out.udb"));
    std::fs::write(&input, basic()).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_mkudb"))
        .args([
            "--word-len",
            "3",
            "-i",
            input.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read(&out).unwrap();
    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&out);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(written, unhex(BASIC_UDB));
}

/// Runs each (binary, arguments) in turn, each on the output of the one
/// before. Returns the last output, which must succeed.
fn pipeline(steps: &[(&str, &[&str])], input: &str) -> Vec<u8> {
    let mut data = input.as_bytes().to_vec();
    for (bin, args) in steps {
        let output = run(bin, args, &data);
        assert!(
            output.status.success(),
            "{bin} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        data = output.stdout;
    }
    data
}

/// Input for the pipelines: lower case, a sequence shorter than 32, a
/// tandem repeat DUST masks, a sequence longer than 50000, and Ns.
fn pipeline_input() -> String {
    fasta(&[
        (
            "p1 x;size=3",
            "GAAGTTGCCGTACTAaattatgacagccggggatcTTCCCGCAAATAGGG",
        ),
        ("p2", "AGGGTCGCAATCGCATCTAA"),
        ("p3", "TTACCACATAGATTCAAGTCATATATATATATATATATATTGCAACCGAT"),
        ("p4", &"ACGT".repeat(12501)),
        ("p5", "CATGTCTACGNNTTGAGAACGTCCNAGACTTGAGTACT"),
    ])
}

const SQFILT: &str = env!("CARGO_BIN_EXE_sqfilt");
const LCMASK: &str = env!("CARGO_BIN_EXE_lcmask");
const MKUDB: &str = env!("CARGO_BIN_EXE_mkudb");
const LENGTHS: &[&str] = &["--min-len", "32", "--max-len", "50000"];

#[test]
fn pipeline_makeudb_usearch() {
    // vsearch --makeudb_usearch in --output out --notrunclabels: DUST
    // masking and the default length filters
    let written = pipeline(
        &[(SQFILT, LENGTHS), (LCMASK, &[]), (MKUDB, &[])],
        &pipeline_input(),
    );
    assert_eq!(written.len(), 262852);
    assert_eq!(sha1(&written), "0a93392b41fbf8c51e97baba8ca308b24711a79d");
}

#[test]
fn pipeline_makeudb_usearch_dbmask_none() {
    // --dbmask none, and --dbmask soft, the same as vsearch upper-cases
    let written = pipeline(
        &[(SQFILT, LENGTHS), (LCMASK, &["--clear"]), (MKUDB, &[])],
        &pipeline_input(),
    );
    assert_eq!(written.len(), 262916);
    assert_eq!(sha1(&written), "5a008e3ccc5f3ff951ee1e895abdb748f75c743f");
}

#[test]
fn pipeline_makeudb_usearch_hardmask() {
    // --dbmask dust --hardmask: upper case, DUST, then masked residues to N
    let written = pipeline(
        &[
            (SQFILT, LENGTHS),
            (LCMASK, &["--clear"]),
            (LCMASK, &["--hard"]),
            (MKUDB, &[]),
        ],
        &pipeline_input(),
    );
    assert_eq!(written.len(), 262852);
    assert_eq!(sha1(&written), "f43d7498bc770c3a73ecd9ec6d95a601217511eb");
}
