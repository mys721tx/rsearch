//! End-to-end tests of the lcatax binary.
//!
//! Expected outputs are those vsearch `--usearch_global --notrunclabels`
//! writes to `--lcaout` beside its `--uc --uc_allhits` hits, from those
//! hits and the database. For each query, it is the taxonomy its hits agree
//! on, level by level, while at least `--cutoff` (`--lca_cutoff`, default
//! 1) of them agree, as vsearch's majority vote picks the candidate name.

use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-lcatax-{}-{name}", std::process::id()))
}

/// Taxonomies differing at every level, in either case, with an empty
/// name, and missing.
const DB: &str = r##">r0;tax=d:B,p:P1,c:C1,g:G1,s:S1
TGGGTAGTGTCACTGCACACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCATTAACTAGCACCAGATGTGATGCGGTTTCATTCCCAGGCCAACAGCAGGACTTGGTCTGAG
>r1;tax=d:B,p:P1,c:C1,g:G1,s:S2
TGGCTAGTGTCACTGCGCACAGTAAACATTCTCGCACATTTTTAACGGGTGAGAGCGCATTAACTATCACCAGATTTGATGCGGTTTCCTGCCCAGGGCAATAGCAGGACTTGGTATGAG
>r2;tax=d:B,p:P1,c:C2,g:G2
TGCCTAGTGTCACTGCGCACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCATTAACTATCACCATATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAGGACTTGGTCTGAG
>r3;tax=d:B,p:P2
TGGCTAGTGTCACTGCGCACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCATTAACTATCACCAGATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAGGACTTGGACTGAG
>r4;tax=D:B,P:P1,c:C1
TGGCTAGTGTCACTGCGCACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCATTAACTATCTCCAGATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAGGACTTGGTCTGAG
>r5;tax=d:B,p:P1,c:,g:G1
TGGCTAGTGCCACTGCGCACAGTAGACATTATCGGACATTTTTAACGGGTGAGCGGGCATTAACTATCACCAGATGTGATGCGGTTTCCTACCCAGGCCAACAGCAGGACTTGGTCTGAG
>r6;tax=d:A,p:P1,c:C1,g:G1,s:S1
TGGCTAGTGTCACTGCGCGCAGTAAACATTATCGCACATGTGTAACTGGTGAGCGGGCATTAACTATCACCAGATGTGATGCGGTTTACTGCACAGGCCAACAGCAGGACTTGGTCTGAG
>r7
TGGCTAGTGTCACTGCGCACAGTAAAAATTATCGCACATTTTTAACGGGTGAGCGGGCATTATCTATCACCAGATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAGGACTTGGTCTGAG
>r8;tax=d:B,p:P1,c:C1,g:G1,s:S1
TGGCTAGTGTCACTGCGCACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCGTTAACTATCACCAGATGTGATGCGGTATCCTGCCCAGGCCAACAGCAGGACTTGCTCTGAG
>r9;tax=d:B,p:P1,c:C1,g:G1,s:S2
TGGCTAGTGTCACTGCGCACAGTAAACAATATCGCGCTTTTTTAACAGGTGAGCGGGCATTAACTACCACCAGATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAGGACTTGGTCTGAG
>r10;tax=d:B,p:P1,c:C2,g:G2
TGGCTAGTGTCACTGCGCACGGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCCTTAACTATCACCAGATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAGGACTTGGTCTGAG
>r11;tax=d:B,p:P2
TGGCTAGTGTCACAGCGCACAGTAAACATTATCGCACATTTATAGCGGGTGACCGGGCATTAACTATCACCAGATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAAGACTTGGTCTGAG
>r12;tax=D:B,P:P1,c:C1
TGGCTTGTGTCACTGCGCACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCATTCACTATCACCAGATGTGATGCGGTTTCCTGCCCAGGCCAACTGCAGGACTTGGTCTGAG
>r13;tax=d:B,p:P1,c:,g:G1
TGGCTAGTGTCACTGCGCGCAGTAAACATTATCCCACATTTTTAACGGGGGATCGGGCATTAACTATCACCAGATGTGATGCGGTTTCCTGCCCAGGCCAACAGCAGAACTTGGTCTGAG
>r14;tax=d:A,p:P1,c:C1,g:G1,s:S1
TGGCTAGTGTCACTGCGCACACTAAACATTATCGCACATTTTTGACGGGTGTGCGGGCATTAATTATCACCAGATGTGATGCGGTTTCCTGCCCAGGCTAACAGCACGACTTGGTCTGAG
>r15
TGGCTAGTGTCACTGCGCACAGTAAACATTATCGCACATTTTTAACGGGTGAGCGGGCATTAACTATCACCAGATGTGATGCGGTTTCCTGCCCGGGCGAACAGCAGGACTAGGTCTGAG
"##;
const HITS: &str = r##"H	2	110	95.5	+	0	0	5I110M5I	q0	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	95.5	+	0	0	5I110M5I	q0	r3;tax=d:B,p:P2
H	4	110	95.5	+	0	0	5I110M5I	q0	r4;tax=D:B,P:P1,c:C1
H	12	110	93.6	+	0	0	5I110M5I	q0	r12;tax=D:B,P:P1,c:C1
H	15	110	93.6	+	0	0	5I110M5I	q0	r15
H	1	110	92.7	+	0	0	5I110M5I	q0	r1;tax=d:B,p:P1,c:C1,g:G1,s:S2
H	2	110	96.4	+	0	0	5I110M5I	q1	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	96.4	+	0	0	5I110M5I	q1	r3;tax=d:B,p:P2
H	4	110	96.4	+	0	0	5I110M5I	q1	r4;tax=D:B,P:P1,c:C1
H	7	110	95.5	+	0	0	5I110M5I	q1	r7
H	10	110	95.5	+	0	0	5I110M5I	q1	r10;tax=d:B,p:P1,c:C2,g:G2
H	8	110	94.5	+	0	0	5I110M5I	q1	r8;tax=d:B,p:P1,c:C1,g:G1,s:S1
H	2	110	95.5	+	0	0	5I110M5I	q2	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	94.5	+	0	0	5I110M5I	q2	r3;tax=d:B,p:P2
H	4	110	94.5	+	0	0	5I110M5I	q2	r4;tax=D:B,P:P1,c:C1
H	7	110	94.5	+	0	0	5I110M5I	q2	r7
H	8	110	92.7	+	0	0	5I110M5I	q2	r8;tax=d:B,p:P1,c:C1,g:G1,s:S1
H	12	110	92.7	+	0	0	5I110M5I	q2	r12;tax=D:B,P:P1,c:C1
H	2	110	97.3	+	0	0	5I110M5I	q3	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	97.3	+	0	0	5I110M5I	q3	r3;tax=d:B,p:P2
H	4	110	97.3	+	0	0	5I110M5I	q3	r4;tax=D:B,P:P1,c:C1
H	7	110	96.4	+	0	0	5I110M5I	q3	r7
H	8	110	95.5	+	0	0	5I110M5I	q3	r8;tax=d:B,p:P1,c:C1,g:G1,s:S1
H	12	110	95.5	+	0	0	5I110M5I	q3	r12;tax=D:B,P:P1,c:C1
H	2	110	95.5	+	0	0	5I110M5I	q4	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	95.5	+	0	0	5I110M5I	q4	r3;tax=d:B,p:P2
H	4	110	95.5	+	0	0	5I110M5I	q4	r4;tax=D:B,P:P1,c:C1
H	7	110	94.5	+	0	0	5I110M5I	q4	r7
H	12	110	93.6	+	0	0	5I110M5I	q4	r12;tax=D:B,P:P1,c:C1
H	15	110	93.6	+	0	0	5I110M5I	q4	r15
H	2	110	99.1	+	0	0	5I110M5I	q5	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q5	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q5	r4;tax=D:B,P:P1,c:C1
H	7	110	98.2	+	0	0	5I110M5I	q5	r7
H	10	110	98.2	+	0	0	5I110M5I	q5	r10;tax=d:B,p:P1,c:C2,g:G2
H	15	110	97.3	+	0	0	5I110M5I	q5	r15
H	2	110	93.6	+	0	0	5I110M5I	q6	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	93.6	+	0	0	5I110M5I	q6	r3;tax=d:B,p:P2
H	4	110	93.6	+	0	0	5I110M5I	q6	r4;tax=D:B,P:P1,c:C1
H	1	110	91.8	+	0	0	5I110M5I	q6	r1;tax=d:B,p:P1,c:C1,g:G1,s:S2
H	12	110	91.8	+	0	0	5I110M5I	q6	r12;tax=D:B,P:P1,c:C1
H	15	110	91.8	+	0	0	5I110M5I	q6	r15
H	2	110	99.1	+	0	0	5I110M5I	q7	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q7	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q7	r4;tax=D:B,P:P1,c:C1
H	7	110	98.2	+	0	0	5I110M5I	q7	r7
H	10	110	98.2	+	0	0	5I110M5I	q7	r10;tax=d:B,p:P1,c:C2,g:G2
H	15	110	97.3	+	0	0	5I110M5I	q7	r15
H	2	110	97.3	+	0	0	5I110M5I	q8	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	97.3	+	0	0	5I110M5I	q8	r3;tax=d:B,p:P2
H	4	110	97.3	+	0	0	5I110M5I	q8	r4;tax=D:B,P:P1,c:C1
H	7	110	96.4	+	0	0	5I110M5I	q8	r7
H	8	110	95.5	+	0	0	5I110M5I	q8	r8;tax=d:B,p:P1,c:C1,g:G1,s:S1
H	15	110	95.5	+	0	0	5I110M5I	q8	r15
H	2	110	99.1	+	0	0	5I110M5I	q9	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q9	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q9	r4;tax=D:B,P:P1,c:C1
H	7	110	98.2	+	0	0	5I110M5I	q9	r7
H	10	110	98.2	+	0	0	5I110M5I	q9	r10;tax=d:B,p:P1,c:C2,g:G2
H	15	110	97.3	+	0	0	5I110M5I	q9	r15
H	2	110	99.1	+	0	0	5I110M5I	q10	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q10	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q10	r4;tax=D:B,P:P1,c:C1
H	7	110	98.2	+	0	0	5I110M5I	q10	r7
H	10	110	98.2	+	0	0	5I110M5I	q10	r10;tax=d:B,p:P1,c:C2,g:G2
H	15	110	97.3	+	0	0	5I110M5I	q10	r15
H	2	110	97.3	+	0	0	5I110M5I	q11	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	97.3	+	0	0	5I110M5I	q11	r3;tax=d:B,p:P2
H	4	110	97.3	+	0	0	5I110M5I	q11	r4;tax=D:B,P:P1,c:C1
H	7	110	96.4	+	0	0	5I110M5I	q11	r7
H	10	110	96.4	+	0	0	5I110M5I	q11	r10;tax=d:B,p:P1,c:C2,g:G2
H	8	110	95.5	+	0	0	5I110M5I	q11	r8;tax=d:B,p:P1,c:C1,g:G1,s:S1
N	*	*	*	.	*	*	*	none	*
"##;

/// Runs lcatax with `args` on `hits` and the database.
fn lcatax(name: &str, hits: &str, args: &[&str]) -> Output {
    let paths = ["hits", "db"].map(|p| temp_path(&format!("{name}-{p}")));
    std::fs::write(&paths[0], hits).unwrap();
    std::fs::write(&paths[1], DB).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_lcatax"))
        .arg(&paths[0])
        .arg("--db")
        .arg(&paths[1])
        .args(args)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    for path in &paths {
        let _ = std::fs::remove_file(path);
    }
    output
}

/// Runs lcatax, which must succeed, and returns what it wrote.
fn lca(name: &str, hits: &str, args: &[&str]) -> String {
    let output = lcatax(name, hits, args);
    assert!(
        output.status.success(),
        "lcatax failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn needs_every_hit_by_default() {
    assert_eq!(
        lca("all", HITS, &[]),
        r##"q0	
q1	
q2	
q3	
q4	
q5	
q6	
q7	
q8	
q9	
q10	
q11	
none	
"##
    );
}

#[test]
fn takes_the_cutoff() {
    assert_eq!(
        lca("c06", HITS, &["--cutoff", "0.6"]),
        r##"q0	d:B,p:P1
q1	d:B,p:P1
q2	d:B,p:P1
q3	d:B,p:P1
q4	d:B
q5	d:B
q6	d:B,p:P1
q7	d:B
q8	d:B
q9	d:B
q10	d:B
q11	d:B,p:P1
none	
"##
    );
    assert_eq!(
        lca("c08", HITS, &["--cutoff", "0.8"]),
        r##"q0	d:B
q1	d:B
q2	d:B
q3	d:B
q4	
q5	
q6	d:B
q7	
q8	
q9	
q10	
q11	d:B
none	
"##
    );
}

#[test]
fn takes_the_top_hits() {
    // the hits of a search with --top_hits_only
    let hits = r##"H	2	110	95.5	+	0	0	5I110M5I	q0	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	95.5	+	0	0	5I110M5I	q0	r3;tax=d:B,p:P2
H	4	110	95.5	+	0	0	5I110M5I	q0	r4;tax=D:B,P:P1,c:C1
H	2	110	96.4	+	0	0	5I110M5I	q1	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	96.4	+	0	0	5I110M5I	q1	r3;tax=d:B,p:P2
H	4	110	96.4	+	0	0	5I110M5I	q1	r4;tax=D:B,P:P1,c:C1
H	2	110	95.5	+	0	0	5I110M5I	q2	r2;tax=d:B,p:P1,c:C2,g:G2
H	2	110	97.3	+	0	0	5I110M5I	q3	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	97.3	+	0	0	5I110M5I	q3	r3;tax=d:B,p:P2
H	4	110	97.3	+	0	0	5I110M5I	q3	r4;tax=D:B,P:P1,c:C1
H	2	110	95.5	+	0	0	5I110M5I	q4	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	95.5	+	0	0	5I110M5I	q4	r3;tax=d:B,p:P2
H	4	110	95.5	+	0	0	5I110M5I	q4	r4;tax=D:B,P:P1,c:C1
H	2	110	99.1	+	0	0	5I110M5I	q5	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q5	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q5	r4;tax=D:B,P:P1,c:C1
H	2	110	93.6	+	0	0	5I110M5I	q6	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	93.6	+	0	0	5I110M5I	q6	r3;tax=d:B,p:P2
H	4	110	93.6	+	0	0	5I110M5I	q6	r4;tax=D:B,P:P1,c:C1
H	2	110	99.1	+	0	0	5I110M5I	q7	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q7	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q7	r4;tax=D:B,P:P1,c:C1
H	2	110	97.3	+	0	0	5I110M5I	q8	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	97.3	+	0	0	5I110M5I	q8	r3;tax=d:B,p:P2
H	4	110	97.3	+	0	0	5I110M5I	q8	r4;tax=D:B,P:P1,c:C1
H	2	110	99.1	+	0	0	5I110M5I	q9	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q9	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q9	r4;tax=D:B,P:P1,c:C1
H	2	110	99.1	+	0	0	5I110M5I	q10	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	99.1	+	0	0	5I110M5I	q10	r3;tax=d:B,p:P2
H	4	110	99.1	+	0	0	5I110M5I	q10	r4;tax=D:B,P:P1,c:C1
H	2	110	97.3	+	0	0	5I110M5I	q11	r2;tax=d:B,p:P1,c:C2,g:G2
H	3	110	97.3	+	0	0	5I110M5I	q11	r3;tax=d:B,p:P2
H	4	110	97.3	+	0	0	5I110M5I	q11	r4;tax=D:B,P:P1,c:C1
N	*	*	*	.	*	*	*	none	*
"##;
    assert_eq!(
        lca("top", hits, &[]),
        r##"q0	d:B
q1	d:B
q2	d:B,p:P1,c:C2,g:G2
q3	d:B
q4	d:B
q5	d:B
q6	d:B
q7	d:B
q8	d:B
q9	d:B
q10	d:B
q11	d:B
none	
"##
    );
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--cutoff", "0.5"][..],
        &["--cutoff", "1.1"],
        &["--cutoff", "x"],
    ] {
        assert_eq!(
            lcatax("args", HITS, args).status.code(),
            Some(2),
            "{args:?}"
        );
    }
    let output = Command::new(env!("CARGO_BIN_EXE_lcatax"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn rejects_hits_not_matching_the_database() {
    let hits = HITS.replacen("H\t", "H\t99", 1);
    assert_eq!(lcatax("far", &hits, &[]).status.code(), Some(1));
}
