//! End-to-end tests of the exsrch binary.
//!
//! Expected hits are those of `vsearch --search_exact --uc --uc_allhits
//! --notrunclabels --minseqlength 0` on the same files. They are every
//! database sequence identical to a query, ignoring case and taking U for
//! T, by target number, the plus strand's first. A query without hits gets
//! an N record. `--maxhits`, `--self`, `--strand` and the abundance filters
//! `--max-qsize`, `--min-tsize`, `--min-size-ratio` and `--max-size-ratio`
//! are vsearch's `--maxhits`, `--self`, `--strand`, `--maxqsize`,
//! `--mintsize`, `--minsizeratio` and `--maxsizeratio`; vsearch's plain
//! `--uc` is `--maxhits 1`. The database is numbered as given. vsearch's
//! `--minseqlength 1`, which drops empty database sequences, is sqfilt.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-exsrch-{}-{name}", std::process::id()))
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

/// Runs exsrch with `args` on `queries` against the database `db`.
fn exsrch(name: &str, db: &str, args: &[&str], queries: &str) -> Output {
    let path = temp_path(&format!("{name}-db"));
    std::fs::write(&path, db).unwrap();
    let mut all = vec!["--db", path.to_str().unwrap()];
    all.extend(args);
    let output = run(env!("CARGO_BIN_EXE_exsrch"), &all, queries.as_bytes());
    let _ = std::fs::remove_file(&path);
    output
}

/// Runs exsrch, which must succeed, and returns its hits.
fn hits(name: &str, db: &str, args: &[&str], queries: &str) -> String {
    let output = exsrch(name, db, args, queries);
    assert!(
        output.status.success(),
        "exsrch failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

/// Identical sequences in either case, with U, an empty one, the reverse
/// complement of the first, and one with N.
const DB: &str = ">d0\nACGTACGTAA\n>d1 x\nacgtacgtaa\n>d2;size=5\nTTACGTACGT\n>d3\n\n\
                  >d4\nACGUACGUAA\n>d5\nACGTACGTAA\n>d6\nACGTNACGTA\n";

const QUERIES: &str = ">q1;size=2\nACGTACGTAA\n>q2\nGGGG\n>empty\n\n>d0\nACGTACGTAA\n\
                       >n\nACGTNACGTA\n>low\nacgtacgtaa\n";

#[test]
fn finds_every_identical_sequence() {
    // an empty query matches an empty sequence
    assert_eq!(
        hits("all", DB, &[], QUERIES),
        r##"H	0	10	100.0	+	0	0	=	q1;size=2	d0
H	1	10	100.0	+	0	0	=	q1;size=2	d1 x
H	4	10	100.0	+	0	0	=	q1;size=2	d4
H	5	10	100.0	+	0	0	=	q1;size=2	d5
N	*	*	*	.	*	*	*	q2	*
H	3	0	100.0	+	0	0	=	empty	d3
H	0	10	100.0	+	0	0	=	d0	d0
H	1	10	100.0	+	0	0	=	d0	d1 x
H	4	10	100.0	+	0	0	=	d0	d4
H	5	10	100.0	+	0	0	=	d0	d5
H	6	10	100.0	+	0	0	=	n	d6
H	0	10	100.0	+	0	0	=	low	d0
H	1	10	100.0	+	0	0	=	low	d1 x
H	4	10	100.0	+	0	0	=	low	d4
H	5	10	100.0	+	0	0	=	low	d5
"##
    );
}

#[test]
fn searches_both_strands() {
    assert_eq!(
        hits("both", DB, &["--strand", "both"], QUERIES),
        r##"H	0	10	100.0	+	0	0	=	q1;size=2	d0
H	1	10	100.0	+	0	0	=	q1;size=2	d1 x
H	2	10	100.0	-	0	0	=	q1;size=2	d2;size=5
H	4	10	100.0	+	0	0	=	q1;size=2	d4
H	5	10	100.0	+	0	0	=	q1;size=2	d5
N	*	*	*	.	*	*	*	q2	*
H	3	0	100.0	+	0	0	=	empty	d3
H	3	0	100.0	-	0	0	=	empty	d3
H	0	10	100.0	+	0	0	=	d0	d0
H	1	10	100.0	+	0	0	=	d0	d1 x
H	2	10	100.0	-	0	0	=	d0	d2;size=5
H	4	10	100.0	+	0	0	=	d0	d4
H	5	10	100.0	+	0	0	=	d0	d5
H	6	10	100.0	+	0	0	=	n	d6
H	0	10	100.0	+	0	0	=	low	d0
H	1	10	100.0	+	0	0	=	low	d1 x
H	2	10	100.0	-	0	0	=	low	d2;size=5
H	4	10	100.0	+	0	0	=	low	d4
H	5	10	100.0	+	0	0	=	low	d5
"##
    );
    // a palindrome matches on both, the plus strand first
    assert_eq!(
        hits("pal", ">t\nACGT\n", &["--strand", "both"], ">p\nACGT\n"),
        r##"H	0	4	100.0	+	0	0	=	p	t
H	0	4	100.0	-	0	0	=	p	t
"##
    );
}

#[test]
fn keeps_the_first_hits() {
    assert_eq!(
        hits("max2", DB, &["--maxhits", "2"], QUERIES),
        r##"H	0	10	100.0	+	0	0	=	q1;size=2	d0
H	1	10	100.0	+	0	0	=	q1;size=2	d1 x
N	*	*	*	.	*	*	*	q2	*
H	3	0	100.0	+	0	0	=	empty	d3
H	0	10	100.0	+	0	0	=	d0	d0
H	1	10	100.0	+	0	0	=	d0	d1 x
H	6	10	100.0	+	0	0	=	n	d6
H	0	10	100.0	+	0	0	=	low	d0
H	1	10	100.0	+	0	0	=	low	d1 x
"##
    );
    // vsearch's --uc without --uc_allhits
    assert_eq!(
        hits("max1", DB, &["--maxhits", "1"], QUERIES),
        r##"H	0	10	100.0	+	0	0	=	q1;size=2	d0
N	*	*	*	.	*	*	*	q2	*
H	3	0	100.0	+	0	0	=	empty	d3
H	0	10	100.0	+	0	0	=	d0	d0
H	6	10	100.0	+	0	0	=	n	d6
H	0	10	100.0	+	0	0	=	low	d0
"##
    );
}

#[test]
fn skips_targets_of_the_same_label() {
    assert_eq!(
        hits("self", DB, &["--self"], QUERIES),
        r##"H	0	10	100.0	+	0	0	=	q1;size=2	d0
H	1	10	100.0	+	0	0	=	q1;size=2	d1 x
H	4	10	100.0	+	0	0	=	q1;size=2	d4
H	5	10	100.0	+	0	0	=	q1;size=2	d5
N	*	*	*	.	*	*	*	q2	*
H	3	0	100.0	+	0	0	=	empty	d3
H	1	10	100.0	+	0	0	=	d0	d1 x
H	4	10	100.0	+	0	0	=	d0	d4
H	5	10	100.0	+	0	0	=	d0	d5
H	6	10	100.0	+	0	0	=	n	d6
H	0	10	100.0	+	0	0	=	low	d0
H	1	10	100.0	+	0	0	=	low	d1 x
H	4	10	100.0	+	0	0	=	low	d4
H	5	10	100.0	+	0	0	=	low	d5
"##
    );
}

/// Queries of abundance 3 and 1 against targets of 4, 1 and 9.
const SIZE_DB: &str = ">a;size=4\nACGTACGTAA\n>b\nACGTACGTAA\n>c;size=9\nACGTACGTAA\n";
const SIZE_QUERIES: &str = ">q;size=3\nACGTACGTAA\n>r\nACGTACGTAA\n";

#[test]
fn filters_by_abundance() {
    let cases = [
        (
            ["--min-tsize", "4"],
            r##"H	0	10	100.0	+	0	0	=	q;size=3	a;size=4
H	2	10	100.0	+	0	0	=	q;size=3	c;size=9
H	0	10	100.0	+	0	0	=	r	a;size=4
H	2	10	100.0	+	0	0	=	r	c;size=9
"##,
        ),
        (
            ["--max-qsize", "2"],
            r##"N	*	*	*	.	*	*	*	q;size=3	*
H	0	10	100.0	+	0	0	=	r	a;size=4
H	1	10	100.0	+	0	0	=	r	b
H	2	10	100.0	+	0	0	=	r	c;size=9
"##,
        ),
        (
            ["--min-size-ratio", "0.5"],
            r##"H	0	10	100.0	+	0	0	=	q;size=3	a;size=4
H	1	10	100.0	+	0	0	=	q;size=3	b
H	1	10	100.0	+	0	0	=	r	b
"##,
        ),
        (
            ["--max-size-ratio", "0.5"],
            r##"H	2	10	100.0	+	0	0	=	q;size=3	c;size=9
H	0	10	100.0	+	0	0	=	r	a;size=4
H	2	10	100.0	+	0	0	=	r	c;size=9
"##,
        ),
    ];
    for (args, expected) in cases {
        assert_eq!(
            hits("size", SIZE_DB, &args, SIZE_QUERIES),
            expected,
            "{args:?}"
        );
    }
}

#[test]
fn reads_fastq_queries() {
    let fastq: String = QUERIES
        .split('>')
        .skip(1)
        .map(|record| {
            let (header, seq) = record.split_once('\n').unwrap();
            let seq = seq.trim_end();
            format!("@{header}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
        })
        .collect();
    assert_eq!(
        hits("fastq", DB, &[], &fastq),
        r##"H	0	10	100.0	+	0	0	=	q1;size=2	d0
H	1	10	100.0	+	0	0	=	q1;size=2	d1 x
H	4	10	100.0	+	0	0	=	q1;size=2	d4
H	5	10	100.0	+	0	0	=	q1;size=2	d5
N	*	*	*	.	*	*	*	q2	*
H	3	0	100.0	+	0	0	=	empty	d3
H	0	10	100.0	+	0	0	=	d0	d0
H	1	10	100.0	+	0	0	=	d0	d1 x
H	4	10	100.0	+	0	0	=	d0	d4
H	5	10	100.0	+	0	0	=	d0	d5
H	6	10	100.0	+	0	0	=	n	d6
H	0	10	100.0	+	0	0	=	low	d0
H	1	10	100.0	+	0	0	=	low	d1 x
H	4	10	100.0	+	0	0	=	low	d4
H	5	10	100.0	+	0	0	=	low	d5
"##
    );
}

#[test]
fn checks_abundances() {
    // the hits of the queries before an invalid abundance are written
    let output = exsrch("zero", ">t\nACGT\n", &[], ">a\nACGT\n>z;size=0\nACGT\n");
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        r##"H	0	4	100.0	+	0	0	=	a	t
"##
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("abundance"), "{stderr}");
    // the database is read whole first
    let output = exsrch("zero-db", ">t;size=0\nACGT\n", &[], QUERIES);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--maxhits", "0"][..],
        &["--maxhits", "x"],
        &["--strand", "minus"],
        &["--min-tsize", "-1"],
        &["--max-qsize", "-1"],
    ] {
        let output = exsrch("args", DB, args, QUERIES);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
    }
    let output = run(env!("CARGO_BIN_EXE_exsrch"), &[], b"");
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn reads_and_writes_files() {
    let (db, input, out) = (
        temp_path("f-db.fa"),
        temp_path("f-in.fa"),
        temp_path("f-out.uc"),
    );
    std::fs::write(&db, DB).unwrap();
    std::fs::write(&input, QUERIES).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_exsrch"))
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
    assert_eq!(
        written,
        r##"H	0	10	100.0	+	0	0	=	q1;size=2	d0
H	1	10	100.0	+	0	0	=	q1;size=2	d1 x
H	4	10	100.0	+	0	0	=	q1;size=2	d4
H	5	10	100.0	+	0	0	=	q1;size=2	d5
N	*	*	*	.	*	*	*	q2	*
H	3	0	100.0	+	0	0	=	empty	d3
H	0	10	100.0	+	0	0	=	d0	d0
H	1	10	100.0	+	0	0	=	d0	d1 x
H	4	10	100.0	+	0	0	=	d0	d4
H	5	10	100.0	+	0	0	=	d0	d5
H	6	10	100.0	+	0	0	=	n	d6
H	0	10	100.0	+	0	0	=	low	d0
H	1	10	100.0	+	0	0	=	low	d1 x
H	4	10	100.0	+	0	0	=	low	d4
H	5	10	100.0	+	0	0	=	low	d5
"##
    );
}
