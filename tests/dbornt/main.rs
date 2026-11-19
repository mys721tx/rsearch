//! End-to-end tests of the dbornt binary.
//!
//! Expected outputs are those of `vsearch --orient --notrunclabels
//! --fasta_width 0` on the same files. Oriented queries go to `-o`
//! (`--fastaout`, or `--fastqout` for FASTQ queries), the others to
//! `--not-oriented` (`--notmatched`). The strand and word counts of every
//! query go to `--tabbedout`, and vsearch's summary to `--report`. Lower
//! case is masked, in the queries as vsearch's default `--qmask dust` does,
//! and in the database as `--dbmask soft` does. A database masked by lcmask
//! is vsearch's default `--dbmask dust`. `--word-len` is `--wordlength`
//! (default 12); a UDB database has its own.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-dbornt-{}-{name}", std::process::id()))
}

/// Builds FASTA text from (label, sequence) pairs.
fn fasta(records: &[(&str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq)| format!(">{id}\n{seq}\n"))
        .collect()
}

/// `n` random bases from xorshift64, seeded by `seed`.
fn random(seed: u64, n: usize) -> String {
    let mut x = seed.wrapping_mul(2_654_435_761).wrapping_add(1);
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGT"[(x >> 62) as usize] as char
        })
        .collect()
}

/// The reverse complement of `seq`, keeping case.
fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            other => other,
        })
        .collect()
}

/// Six random sequences of 200 bases, which DUST leaves alone.
fn database() -> Vec<String> {
    (0..6).map(|i| random(i, 200)).collect()
}

fn database_fasta(db: &[String]) -> String {
    db.iter()
        .enumerate()
        .map(|(i, s)| format!(">d{i}\n{s}\n"))
        .collect()
}

/// Queries: forward and reverse parts of the database, one of neither, one
/// in lower case, half and half, an empty one, a reverse one partly in
/// lower case, and one shorter than the words.
fn queries() -> Vec<(String, String)> {
    let db = database();
    let part = rc(&db[5][100..180]);
    [
        ("fwd", db[0][20..120].to_string()),
        ("rev desc", rc(&db[1][50..150])),
        ("none", random(99, 100)),
        ("low", db[2][..80].to_lowercase()),
        ("mixed", format!("{}{}", &db[3][..50], rc(&db[4][..50]))),
        ("empty", String::new()),
        (
            "part;size=3",
            format!("{}{}", part[..30].to_lowercase(), &part[30..]),
        ),
        ("short", "ACGT".to_string()),
    ]
    .into_iter()
    .map(|(h, s)| (h.to_string(), s))
    .collect()
}

fn queries_fasta() -> String {
    queries()
        .iter()
        .map(|(h, s)| format!(">{h}\n{s}\n"))
        .collect()
}

/// Scores that differ along a read, so that reversing them shows.
fn queries_fastq() -> String {
    queries()
        .iter()
        .map(|(h, s)| {
            let qual: String = (0..s.len())
                .map(|i| (33 + (i % 40) as u8) as char)
                .collect();
            format!("@{h}\n{s}\n+\n{qual}\n")
        })
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
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&input);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// What a run wrote.
struct Run {
    output: Output,
    oriented: String,
    not_oriented: String,
    table: String,
    report: String,
}

/// Runs dbornt with `args` on `input`, fed on stdin, against the database
/// `db` (the bytes of a FASTA, FASTQ or UDB file), writing every output.
/// Temporary files are named after `name`.
fn dbornt(name: &str, db: &[u8], args: &[&str], input: &str) -> Run {
    let paths = ["db", "not", "table", "report"].map(|p| temp_path(&format!("{name}-{p}")));
    std::fs::write(&paths[0], db).unwrap();
    let path = |i: usize| paths[i].to_str().unwrap().to_string();
    let mut all = vec![
        "--db".to_string(),
        path(0),
        "--not-oriented".to_string(),
        path(1),
        "--tabbedout".to_string(),
        path(2),
        "--report".to_string(),
        path(3),
    ];
    all.extend(args.iter().map(|a| a.to_string()));
    let all: Vec<&str> = all.iter().map(String::as_str).collect();
    let output = run(env!("CARGO_BIN_EXE_dbornt"), &all, input.as_bytes());
    let read = |i: usize| std::fs::read_to_string(&paths[i]).unwrap_or_default();
    let result = Run {
        oriented: String::from_utf8(output.stdout.clone()).unwrap(),
        not_oriented: read(1),
        table: read(2),
        report: read(3),
        output,
    };
    for path in &paths {
        let _ = std::fs::remove_file(path);
    }
    result
}

/// Runs dbornt, which must succeed.
fn dbornt_ok(name: &str, db: &[u8], args: &[&str], input: &str) -> Run {
    let result = dbornt(name, db, args, input);
    assert!(
        result.output.status.success(),
        "dbornt failed: {}",
        String::from_utf8_lossy(&result.output.stderr)
    );
    result
}

const TABLE: &str = r##"fwd	+	89	0
rev desc	-	0	89
none	?	0	0
low	?	0	0
mixed	?	39	39
empty	?	0	0
part;size=3	-	0	39
short	?	0	0
"##;

const ORIENTED: &str = r##">fwd
TCACTTCACCGTAAGGATGAATACTTAATGACACGCTATCGCAAGACGTCACAGTTTATGACAGACGGACCTAGTTTACGCAACCAAACATCTAGTACGG
>rev desc
TGTATGTATTAATCTGTGTTCAAGGGCTTCGATACTCAATAAGGAACCAAGATCCAGGCTCATAATCCACCCTCCCAATAGGGAAGTAGGGCGTCAATCG
>part;size=3
GTAGGAAAGCCCAGCAAGCAGATCATGGCACAATTCATAGGTATTGGAAGacccctcagcagcccccgtttagtattatc
"##;

#[test]
fn orients_queries() {
    let db = database_fasta(&database());
    let result = dbornt_ok("fasta", db.as_bytes(), &[], &queries_fasta());
    assert_eq!(result.oriented, ORIENTED);
    assert_eq!(
        result.not_oriented,
        r##">none
CCTCTCTGCGGTAAGGCGCTGAGAGCCAGGATCAATCGACCGCTCATACCAATTCTTCTTTCATCCGGGTCAAAGAACTGTATAGTGTGCTATCTGTCTC
>low
cgtctcggaatactctcctactgcgccagcgtgcttccgaaggcagtgacgtcggccgtccccacaactagctccgtgca
>mixed
CAGATGACAGCGCTCTTAAGCCCCACGCATGTGTCAGTAGGCCCCACCGTGTTCTGACACTGCCGGAAGACCCTGAATAGATCTCGGTTTCGGTCGGATC
>empty

>short
ACGT
"##
    );
    assert_eq!(result.table, TABLE);
    assert_eq!(
        result.report,
        "Forward oriented sequences: 1 (12.50%)\nReverse oriented sequences: 2 \
         (25.00%)\nAll oriented sequences:     3 (37.50%)\nNot oriented sequences:     \
         5 (62.50%)\nTotal number of sequences:  8\n"
    );
    assert!(result.output.stderr.is_empty());
}

#[test]
fn orients_fastq_queries() {
    // the scores of a reversed read are reversed
    let db = database_fasta(&database());
    let result = dbornt_ok("fastq", db.as_bytes(), &[], &queries_fastq());
    assert_eq!(
        result.oriented,
        r##"@fwd
TCACTTCACCGTAAGGATGAATACTTAATGACACGCTATCGCAAGACGTCACAGTTTATGACAGACGGACCTAGTTTACGCAACCAAACATCTAGTACGG
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH!"#$%&'()*+,-./01234
@rev desc
TGTATGTATTAATCTGTGTTCAAGGGCTTCGATACTCAATAAGGAACCAAGATCCAGGCTCATAATCCACCCTCCCAATAGGGAAGTAGGGCGTCAATCG
+
43210/.-,+*)('&%$#"!HGFEDCBA@?>=<;:9876543210/.-,+*)('&%$#"!HGFEDCBA@?>=<;:9876543210/.-,+*)('&%$#"!
@part;size=3
GTAGGAAAGCCCAGCAAGCAGATCATGGCACAATTCATAGGTATTGGAAGacccctcagcagcccccgtttagtattatc
+
HGFEDCBA@?>=<;:9876543210/.-,+*)('&%$#"!HGFEDCBA@?>=<;:9876543210/.-,+*)('&%$#"!
"##
    );
    assert_eq!(
        result.not_oriented,
        r##"@none
CCTCTCTGCGGTAAGGCGCTGAGAGCCAGGATCAATCGACCGCTCATACCAATTCTTCTTTCATCCGGGTCAAAGAACTGTATAGTGTGCTATCTGTCTC
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH!"#$%&'()*+,-./01234
@low
cgtctcggaatactctcctactgcgccagcgtgcttccgaaggcagtgacgtcggccgtccccacaactagctccgtgca
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH
@mixed
CAGATGACAGCGCTCTTAAGCCCCACGCATGTGTCAGTAGGCCCCACCGTGTTCTGACACTGCCGGAAGACCCTGAATAGATCTCGGTTTCGGTCGGATC
+
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGH!"#$%&'()*+,-./01234
@empty

+

@short
ACGT
+
!"#$
"##
    );
    assert_eq!(result.table, TABLE);
}

const TABLE_8: &str = r##"fwd	+	91	0
rev desc	-	0	91
none	-	0	4
low	?	0	0
mixed	?	43	41
empty	?	0	0
part;size=3	-	0	43
short	?	0	0
"##;

#[test]
fn takes_the_word_length() {
    let db = database_fasta(&database());
    let result = dbornt_ok(
        "len8",
        db.as_bytes(),
        &["--word-len", "8"],
        &queries_fasta(),
    );
    assert_eq!(result.table, TABLE_8);
}

/// The UDB file mkudb writes for `fasta`, at word length 8.
fn udb(fasta: &str) -> Vec<u8> {
    let output = run(env!("CARGO_BIN_EXE_mkudb"), &[], fasta.as_bytes());
    assert!(output.status.success());
    output.stdout
}

#[test]
fn reads_udb_databases() {
    // the word length is the file's, with vsearch's warning when it is not
    // the one asked for
    let db = udb(&database_fasta(&database()));
    let result = dbornt_ok("udb", &db, &[], &queries_fasta());
    assert_eq!(
        result.table,
        r##"fwd	+	91	0
rev desc	-	0	91
none	-	0	4
low	?	0	0
mixed	?	43	41
empty	?	0	0
part;size=3	-	0	43
short	?	0	0
"##
    );
    assert_eq!(result.table, TABLE_8);
    assert_eq!(
        String::from_utf8_lossy(&result.output.stderr),
        "Warning: the word length of the UDB file, 8, is used\n"
    );
    let result = dbornt_ok("udb8", &db, &["--word-len", "8"], &queries_fasta());
    assert_eq!(result.table, TABLE_8);
    assert!(result.output.stderr.is_empty());
}

#[test]
fn skips_lowercase_database_words() {
    // vsearch --dbmask soft: the first sequence all masked, the second half
    let mut db = database();
    db[0] = db[0].to_lowercase();
    db[1] = format!("{}{}", db[1][..100].to_lowercase(), &db[1][100..]);
    let result = dbornt_ok(
        "lower",
        database_fasta(&db).as_bytes(),
        &[],
        &queries_fasta(),
    );
    assert_eq!(
        result.table,
        r##"fwd	?	0	0
rev desc	-	0	39
none	?	0	0
low	?	0	0
mixed	?	39	39
empty	?	0	0
part;size=3	-	0	39
short	?	0	0
"##
    );
    let stderr = String::from_utf8_lossy(&result.output.stderr);
    assert!(
        stderr.starts_with("Warning: 1 of 6 sequences give no word"),
        "{stderr}"
    );
}

/// The table line of query `AAC` against a database of `fwd` sequences
/// `AAC` and `rev` sequences `GTT`, at word length 3.
fn strand(fwd: usize, rev: usize) -> String {
    let mut db = String::new();
    for i in 0..fwd {
        db.push_str(&format!(">f{i}\nAAC\n"));
    }
    for i in 0..rev {
        db.push_str(&format!(">r{i}\nGTT\n"));
    }
    dbornt_ok("strand", db.as_bytes(), &["--word-len", "3"], ">q\nAAC\n").table
}

#[test]
fn needs_eight_times_the_hits_on_one_strand() {
    assert_eq!(strand(9, 1), "q\t+\t1\t0\n");
    assert_eq!(strand(8, 1), "q\t?\t0\t0\n");
    assert_eq!(strand(1, 0), "q\t+\t1\t0\n");
    assert_eq!(strand(0, 1), "q\t-\t0\t1\n");
    assert_eq!(strand(0, 0), "q\t?\t0\t0\n");
    assert_eq!(strand(17, 2), "q\t+\t1\t0\n");
    assert_eq!(strand(16, 2), "q\t?\t0\t0\n");
}

#[test]
fn needs_four_times_the_words_on_one_strand() {
    // AAC, AAG, ACA, AGA and CAA are forward words; ATG, ATC, ACC and AGC
    // reverse ones; each word counts once
    let db = fasta(&[
        ("a", "AAC"),
        ("b", "AAG"),
        ("c", "ACA"),
        ("d", "AGA"),
        ("e", "CAA"),
        ("f", "CAT"),
        ("g", "GAT"),
        ("h", "GGT"),
        ("i", "GCT"),
        ("p", "ACGT"),
    ]);
    let queries = fasta(&[
        ("q41", "AACNAAGNACANAGANATG"),
        ("q31", "AACNAAGNACANATG"),
        ("q14", "AACNATGNATCNACCNAGC"),
        ("q13", "AACNATGNATCNACC"),
        ("q10", "AACNAACNAAC"),
        ("q01", "ATG"),
    ]);
    let result = dbornt_ok("words", db.as_bytes(), &["--word-len", "3"], &queries);
    assert_eq!(
        result.table,
        "q41\t+\t4\t1\nq31\t?\t3\t1\nq14\t-\t1\t4\nq13\t?\t1\t3\nq10\t+\t1\t0\nq01\t-\t0\t1\n"
    );
    // a word that is its own reverse complement counts for neither strand
    let result = dbornt_ok(
        "palindrome",
        db.as_bytes(),
        &["--word-len", "4"],
        ">pal\nACGT\n",
    );
    assert_eq!(result.table, "pal\t?\t0\t0\n");
}

#[test]
fn orients_nothing_against_an_empty_database() {
    let result = dbornt_ok("empty-db", b"", &[], &queries_fasta());
    assert_eq!(
        result.table,
        r##"fwd	?	0	0
rev desc	?	0	0
none	?	0	0
low	?	0	0
mixed	?	0	0
empty	?	0	0
part;size=3	?	0	0
short	?	0	0
"##
    );
    assert_eq!(result.oriented, "");
    assert_eq!(result.not_oriented, queries_fasta());
    assert_eq!(
        result.report,
        "Forward oriented sequences: 0 (0.00%)\nReverse oriented sequences: 0 \
         (0.00%)\nAll oriented sequences:     0 (0.00%)\nNot oriented sequences:     \
         8 (100.00%)\nTotal number of sequences:  8\n"
    );
}

#[test]
fn reports_no_queries() {
    let db = database_fasta(&database());
    let result = dbornt_ok("no-queries", db.as_bytes(), &[], "");
    assert_eq!(
        result.report,
        "Forward oriented sequences: 0\nReverse oriented sequences: 0\nAll oriented \
         sequences:     0\nNot oriented sequences:     0\nTotal number of sequences:  0\n"
    );
    assert_eq!(result.table, "");
}

#[test]
fn checks_abundances() {
    // the database is read whole first
    let db = fasta(&[("d;size=0", "ACGT")]);
    let result = dbornt("zero-db", db.as_bytes(), &[], &queries_fasta());
    assert_eq!(result.output.status.code(), Some(1));
    assert_eq!(result.oriented, "");
    assert_eq!(result.table, "");
    // the queries before an invalid abundance are written
    let db = database_fasta(&database());
    let input = queries_fasta().replace(">rev desc", ">rev desc;size=0");
    let result = dbornt("zero-query", db.as_bytes(), &[], &input);
    assert_eq!(result.output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&result.output.stderr);
    assert!(stderr.contains("abundance"), "{stderr}");
    assert_eq!(result.oriented, format!(">fwd\n{}\n", queries()[0].1));
    assert_eq!(result.table, "fwd\t+\t89\t0\n");
}

#[test]
fn rejects_invalid_input() {
    // the queries' format is checked before the database is read
    let db = fasta(&[("d;size=0", "ACGT")]);
    let result = dbornt("bad-query", db.as_bytes(), &[], "xx\n");
    assert_eq!(result.output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&result.output.stderr),
        "Error: File type not recognized.\n"
    );
    // a truncated UDB database
    let db = udb(&database_fasta(&database()));
    let result = dbornt(
        "bad-udb",
        &db[..300],
        &["--word-len", "8"],
        &queries_fasta(),
    );
    assert_eq!(result.output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&result.output.stderr),
        "Error: invalid UDB file: truncated\n"
    );
    assert_eq!(result.oriented, "");
}

#[test]
fn rejects_invalid_arguments() {
    for args in [&["--word-len", "2"][..], &["--word-len", "16"]] {
        let result = dbornt("args", b">d\nACGT\n", args, ">q\nACGT\n");
        assert_eq!(result.output.status.code(), Some(2), "{args:?}");
    }
    // --db is required
    let output = Command::new(env!("CARGO_BIN_EXE_dbornt"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn reads_and_writes_files() {
    let (db, input, out) = (
        temp_path("f-db.fa"),
        temp_path("f-in.fa"),
        temp_path("f-out.fa"),
    );
    std::fs::write(&db, database_fasta(&database())).unwrap();
    std::fs::write(&input, queries_fasta()).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dbornt"))
        .args(["--db", db.to_str().unwrap(), "-i", input.to_str().unwrap()])
        .args(["-o", out.to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    for path in [&db, &input, &out] {
        let _ = std::fs::remove_file(path);
    }
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(written, ORIENTED);
}

#[test]
fn pipeline_orient() {
    // vsearch --orient q --db db --fastaout out: DUST on the database, which
    // lcmask does first, and FASTA lines of 80, which conv writes
    let db = temp_path("pipe-db.fa");
    let masked = run(
        env!("CARGO_BIN_EXE_lcmask"),
        &[],
        database_fasta(&database()).as_bytes(),
    );
    assert!(masked.status.success());
    std::fs::write(&db, masked.stdout).unwrap();
    let oriented = run(
        env!("CARGO_BIN_EXE_dbornt"),
        &["--db", db.to_str().unwrap()],
        queries_fastq().as_bytes(),
    );
    let _ = std::fs::remove_file(&db);
    assert!(oriented.status.success());
    let written = run(env!("CARGO_BIN_EXE_conv"), &["--fasta"], &oriented.stdout);
    assert_eq!(
        String::from_utf8(written.stdout).unwrap(),
        r##">fwd
TCACTTCACCGTAAGGATGAATACTTAATGACACGCTATCGCAAGACGTCACAGTTTATGACAGACGGACCTAGTTTACG
CAACCAAACATCTAGTACGG
>rev desc
TGTATGTATTAATCTGTGTTCAAGGGCTTCGATACTCAATAAGGAACCAAGATCCAGGCTCATAATCCACCCTCCCAATA
GGGAAGTAGGGCGTCAATCG
>part;size=3
GTAGGAAAGCCCAGCAAGCAGATCATGGCACAATTCATAGGTATTGGAAGacccctcagcagcccccgtttagtattatc
"##
    );
}
