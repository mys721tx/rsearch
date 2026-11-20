//! End-to-end tests of the otutab binary.
//!
//! Expected tables are those vsearch (`--search_exact` or
//! `--usearch_global`, with `--notrunclabels`) writes to `--otutabout`,
//! `--mothur_shared_out` (`--mothur`) and `--biomout` (`--biom`) beside its
//! `--uc --uc_allhits` hits, from those hits and the database. Each query
//! counts its abundance for the sample of its label and the OTU of its
//! first hit's target. Every database sequence is an OTU. In BIOM, the
//! lines naming the file, the program and the date are rsearch's own.

use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-otutab-{}-{name}", std::process::id()))
}

/// Runs otutab with `args` on `hits` and the database `db`.
fn otutab(name: &str, hits: &str, db: &str, args: &[&str]) -> Output {
    let paths = ["hits", "db"].map(|p| temp_path(&format!("{name}-{p}")));
    std::fs::write(&paths[0], hits).unwrap();
    std::fs::write(&paths[1], db).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_otutab"))
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

/// Runs otutab, which must succeed, and returns its table.
fn table(name: &str, hits: &str, db: &str, args: &[&str]) -> String {
    let output = otutab(name, hits, db, args);
    assert!(
        output.status.success(),
        "otutab failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

/// the BIOM table without the lines naming the file, the program and the
/// date.
fn biom(text: &str) -> String {
    text.lines()
        .filter(|l| {
            !l.starts_with("\t\"id\":")
                && !l.contains("\"generated_by\"")
                && !l.contains("\"date\"")
        })
        .map(|l| format!("{l}\n"))
        .collect()
}

/// Samples from sample= and barcodelabel= annotations, whichever comes
/// first, or the label's leading letters, digits and underscores; OTUs from
/// otu= annotations or the label up to its first semicolon; abundances.
const EXACT_DB: &str = r##">t1;otu=A;tax=d:X
ACGTACGTAAGGCCTTACGT
>t2
GGGGCCCCAAAATTTTACGA
>t3;tax=d:Y
TTTTGGGGCCCCAAAATACG
>t4;x
TTTTGGGGCCCCAAAATACG
"##;
const EXACT_HITS: &str = r##"H	0	20	100.0	+	0	0	=	s1.r1;size=5	t1;otu=A;tax=d:X
H	1	20	100.0	+	0	0	=	s1.r2;sample=zz	t2
H	0	20	100.0	+	0	0	=	s2_x;size=2	t1;otu=A;tax=d:X
N	*	*	*	.	*	*	*	s2;barcodelabel=bb;sample=cc	*
"##;

#[test]
fn counts_samples_and_otus() {
    assert_eq!(
        table("exact", EXACT_HITS, EXACT_DB, &[]),
        r##"#OTU ID	bb	s1	s2_x	zz	taxonomy
A	0	5	2	0	d:X
t2	0	0	0	1	
t3	0	0	0	0	d:Y
t4	0	0	0	0	
"##
    );
    assert_eq!(
        table("exact-mothur", EXACT_HITS, EXACT_DB, &["--mothur"]),
        r##"label	Group	numOtus	A	t2	t3	t4
vsearch	bb	4	0	0	0	0
vsearch	s1	4	5	0	0	0
vsearch	s2_x	4	2	0	0	0
vsearch	zz	4	0	1	0	0
"##
    );
    let written = table("exact-biom", EXACT_HITS, EXACT_DB, &["--biom"]);
    assert_eq!(
        biom(&written),
        r##"{
	"format": "Biological Observation Matrix 1.0",
	"format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
	"type": "OTU table",
	"matrix_type": "sparse",
	"matrix_element_type": "int",
	"shape": [4,4],
	"rows":[
		{"id":"A", "metadata":{"taxonomy":"d:X"}},
		{"id":"t2", "metadata":{"taxonomy":""}},
		{"id":"t3", "metadata":{"taxonomy":"d:Y"}},
		{"id":"t4", "metadata":{"taxonomy":""}}
	],
	"columns":[
		{"id":"bb", "metadata":null},
		{"id":"s1", "metadata":null},
		{"id":"s2_x", "metadata":null},
		{"id":"zz", "metadata":null}
	],
	"data": [
		[0,1,5],
		[0,2,2],
		[1,3,1]
	]
}
"##
    );
    assert!(
        written.contains("\t\"generated_by\": \"rsearch "),
        "{written}"
    );
}

/// A search with weak hits and several hits a query, both strands, OTU
/// names shared by two sequences.
const SEARCH_DB: &str = r##">otuA;tax=d:B,p:F
CGATTCAAATGACGGCAGCAGGCCGGGAGTCCCTGAGAGGCTTGTTCCGGAAATGTGCCATCTGCGTGCGAACGCAGCGTAAGAGGAGGGCTAGCTGCGT
>otuB;size=3;tax=d:B,p:P
CGAGATCGGGATCTCAAAACCATCGAAGTCTCCTTTACTTCTCTCAAGGCCCTGCGAGATATTATCCGGTGTCGGTTAGCATCGACTTTTCACCAGATTC
>x;otu=otuC
ACCGTTAAAATGCAGAAGGAATTCGTCTTAAAGTTTACGTTACGCCCGTGGACAGAATTACTGGCCAAGTGTTTCGGGCTACCGGCGAATCGGGCGAAAG
>otuA;tax=d:B,p:G
CGATTCAAATGACGACAGCAGGCCGGGAGTCACCGAGGGGCTTGTTCCGGAAATGTGCCATCTGCGTGCGAACGCAGCGTGAGAAGCGGGCCAGGCACGT
>otuD
ACCTAACTCGTCTCGGCGTTTGGTTAGCTGTTACATGGAGTCTACCAGCACTGACCACAGAAGAGCCTCAGACTCCGTATCCGTGTGCTATAGAGCACAA
"##;
const SEARCH_HITS: &str = r##"N	*	*	*	.	*	*	*	r0;barcodelabel=S3;size=2	*
N	*	*	*	.	*	*	*	S2.1	*
H	1	90	96.7	+	0	0	5I90M5I	S1_2;size=5	otuB;size=3;tax=d:B,p:P
N	*	*	*	.	*	*	*	r3;sample=S1	*
H	4	90	94.4	+	0	0	5I90M5I	r4;sample=S1	otuD
H	1	90	96.7	+	0	0	5I90M5I	r5;sample=S1	otuB;size=3;tax=d:B,p:P
H	2	90	100.0	+	0	0	5I90M5I	S2_6;size=8	x;otu=otuC
H	1	90	81.1	+	0	0	5I90M5I	r7;sample=S1	otuB;size=3;tax=d:B,p:P
H	4	90	98.9	+	0	0	5I90M5I	S2.8	otuD
N	*	*	*	.	*	*	*	r9;sample=S3	*
H	0	90	93.3	+	0	0	5I90M5I	S3_10;size=8	otuA;tax=d:B,p:F
H	3	90	83.3	+	0	0	5I90M5I	S3_10;size=8	otuA;tax=d:B,p:G
H	1	90	93.3	+	0	0	5I90M5I	r11;sample=S2	otuB;size=3;tax=d:B,p:P
H	4	90	100.0	+	0	0	5I90M5I	S2.12	otuD
H	0	90	100.0	+	0	0	5I90M5I	r13;sample=S3	otuA;tax=d:B,p:F
H	3	90	90.0	+	0	0	5I90M5I	r13;sample=S3	otuA;tax=d:B,p:G
H	1	90	97.8	+	0	0	5I90M5I	S2_14;size=7	otuB;size=3;tax=d:B,p:P
N	*	*	*	.	*	*	*	S2.15	*
H	1	90	100.0	+	0	0	5I90M5I	r16;barcodelabel=S2;size=2	otuB;size=3;tax=d:B,p:P
H	1	90	82.2	+	0	0	5I90M5I	S2.17	otuB;size=3;tax=d:B,p:P
H	2	90	85.6	+	0	0	5I90M5I	r18;barcodelabel=S1;size=2	x;otu=otuC
H	4	90	95.6	+	0	0	5I90M5I	S1_19;size=6	otuD
H	4	90	97.8	+	0	0	5I90M5I	r20;sample=S3	otuD
H	4	90	100.0	+	0	0	5I90M5I	r21;sample=S3	otuD
H	1	90	100.0	+	0	0	5I90M5I	r22;sample=S1	otuB;size=3;tax=d:B,p:P
H	1	90	100.0	+	0	0	5I90M5I	S3_23;size=8	otuB;size=3;tax=d:B,p:P
H	2	90	93.3	+	0	0	5I90M5I	r24;barcodelabel=S1;size=2	x;otu=otuC
H	2	90	100.0	+	0	0	5I90M5I	r25;sample=S3	x;otu=otuC
H	4	90	82.2	+	0	0	5I90M5I	S3_26;size=2	otuD
H	0	90	100.0	+	0	0	5I90M5I	S3.27	otuA;tax=d:B,p:F
H	3	90	90.0	+	0	0	5I90M5I	S3.27	otuA;tax=d:B,p:G
N	*	*	*	.	*	*	*	r28;sample=S2	*
H	1	90	100.0	+	0	0	5I90M5I	S2_29;size=7	otuB;size=3;tax=d:B,p:P
"##;

#[test]
fn counts_the_first_hit_of_each_query() {
    assert_eq!(
        table("search", SEARCH_HITS, SEARCH_DB, &[]),
        r##"#OTU ID	S1	S1_19	S1_2	S2	S2_14	S2_29	S2_6	S3	S3_10	S3_23	S3_26	taxonomy
otuA	0	0	0	0	0	0	0	2	8	0	0	d:B,p:F
otuB	3	0	5	4	7	7	0	0	0	8	0	d:B,p:P
otuC	4	0	0	0	0	0	8	1	0	0	0	
otuD	1	6	0	2	0	0	0	2	0	0	2	
"##
    );
    assert_eq!(
        table("search-mothur", SEARCH_HITS, SEARCH_DB, &["--mothur"]),
        r##"label	Group	numOtus	otuA	otuB	otuC	otuD
vsearch	S1	4	0	3	4	1
vsearch	S1_19	4	0	0	0	6
vsearch	S1_2	4	0	5	0	0
vsearch	S2	4	0	4	0	2
vsearch	S2_14	4	0	7	0	0
vsearch	S2_29	4	0	7	0	0
vsearch	S2_6	4	0	0	8	0
vsearch	S3	4	2	0	1	2
vsearch	S3_10	4	8	0	0	0
vsearch	S3_23	4	0	8	0	0
vsearch	S3_26	4	0	0	0	2
"##
    );
    let written = table("search-biom", SEARCH_HITS, SEARCH_DB, &["--biom"]);
    assert_eq!(
        biom(&written),
        r##"{
	"format": "Biological Observation Matrix 1.0",
	"format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
	"type": "OTU table",
	"matrix_type": "sparse",
	"matrix_element_type": "int",
	"shape": [4,11],
	"rows":[
		{"id":"otuA", "metadata":{"taxonomy":"d:B,p:F"}},
		{"id":"otuB", "metadata":{"taxonomy":"d:B,p:P"}},
		{"id":"otuC", "metadata":{"taxonomy":""}},
		{"id":"otuD", "metadata":{"taxonomy":""}}
	],
	"columns":[
		{"id":"S1", "metadata":null},
		{"id":"S1_19", "metadata":null},
		{"id":"S1_2", "metadata":null},
		{"id":"S2", "metadata":null},
		{"id":"S2_14", "metadata":null},
		{"id":"S2_29", "metadata":null},
		{"id":"S2_6", "metadata":null},
		{"id":"S3", "metadata":null},
		{"id":"S3_10", "metadata":null},
		{"id":"S3_23", "metadata":null},
		{"id":"S3_26", "metadata":null}
	],
	"data": [
		[0,7,2],
		[0,8,8],
		[1,0,3],
		[1,2,5],
		[1,3,4],
		[1,4,7],
		[1,5,7],
		[1,9,8],
		[2,0,4],
		[2,6,8],
		[2,7,1],
		[3,0,1],
		[3,1,6],
		[3,3,2],
		[3,7,2],
		[3,10,2]
	]
}
"##
    );
}

#[test]
fn takes_empty_taxonomies() {
    // an empty tax= still makes a taxonomy column
    let db = r##">z1;tax=
ACGT
>z2
ACGG
"##;
    let hits = r##"H	0	4	100.0	+	0	0	=	n1	z1;tax=
"##;
    assert_eq!(
        table("empty-tax", hits, db, &[]),
        r##"#OTU ID	n1	taxonomy
z1	1	
z2	0	
"##
    );
    assert_eq!(
        table("empty-tax-mothur", hits, db, &["--mothur"]),
        r##"label	Group	numOtus	z1	z2
vsearch	n1	2	1	0
"##
    );
    let written = table("empty-tax-biom", hits, db, &["--biom"]);
    assert_eq!(
        biom(&written),
        r##"{
	"format": "Biological Observation Matrix 1.0",
	"format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
	"type": "OTU table",
	"matrix_type": "sparse",
	"matrix_element_type": "int",
	"shape": [2,1],
	"rows":[
		{"id":"z1", "metadata":{"taxonomy":""}},
		{"id":"z2", "metadata":{"taxonomy":""}}
	],
	"columns":[
		{"id":"n1", "metadata":null}
	],
	"data": [
		[0,0,1]
	]
}
"##
    );
}

#[test]
fn lists_every_otu_without_queries() {
    assert_eq!(
        table("none", "", EXACT_DB, &[]),
        r##"#OTU ID	taxonomy
A	d:X
t2	
t3	d:Y
t4	
"##
    );
    assert_eq!(
        table("none-mothur", "", EXACT_DB, &["--mothur"]),
        r##"label	Group	numOtus	A	t2	t3	t4
"##
    );
    let written = table("none-biom", "", EXACT_DB, &["--biom"]);
    assert_eq!(
        biom(&written),
        r##"{
	"format": "Biological Observation Matrix 1.0",
	"format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
	"type": "OTU table",
	"matrix_type": "sparse",
	"matrix_element_type": "int",
	"shape": [4,0],
	"rows":[
		{"id":"A", "metadata":{"taxonomy":"d:X"}},
		{"id":"t2", "metadata":{"taxonomy":""}},
		{"id":"t3", "metadata":{"taxonomy":"d:Y"}},
		{"id":"t4", "metadata":{"taxonomy":""}}
	],
	"columns":[
	],
	"data": [
	]
}
"##
    );
}

#[test]
fn rejects_hits_not_matching_the_database() {
    let hits = EXACT_HITS.replacen("H\t0\t", "H\t9\t", 1);
    assert_eq!(otutab("far", &hits, EXACT_DB, &[]).status.code(), Some(1));
    assert_eq!(
        otutab("bad", "X\tq\n", EXACT_DB, &[]).status.code(),
        Some(1)
    );
}

#[test]
fn rejects_invalid_arguments() {
    for args in [&["--mothur", "--biom"][..]] {
        assert_eq!(
            otutab("args", EXACT_HITS, EXACT_DB, args).status.code(),
            Some(2)
        );
    }
    let output = Command::new(env!("CARGO_BIN_EXE_otutab"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn writes_files() {
    let (hits, db, out) = (temp_path("f-hits"), temp_path("f-db"), temp_path("f-out"));
    std::fs::write(&hits, EXACT_HITS).unwrap();
    std::fs::write(&db, EXACT_DB).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_otutab"))
        .arg(&hits)
        .args([
            "--db",
            db.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--biom",
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    for path in [&hits, &db, &out] {
        let _ = std::fs::remove_file(path);
    }
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    // the BIOM table's id is the output file
    assert!(written.contains(&format!("\t\"id\":\"{}\",", out.to_str().unwrap())));
}
