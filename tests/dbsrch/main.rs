//! End-to-end tests of the dbsrch binary.
//!
//! Expected hits are those of `vsearch --usearch_global --uc --uc_allhits
//! --notrunclabels --qmask soft --dbmask soft --minseqlength 1` on the same
//! queries and database. The targets sharing the most words with a query
//! are aligned, best first, by vsearch's 16-bit aligner (or its
//! linear-memory one), until `--maxaccepts` are accepted or `--maxrejects`
//! rejected. The hits accepted, and the weak ones, are written by
//! decreasing identity, then target. The options are vsearch's, in words:
//! `--word-len` is `--wordlength`, `--min-word-matches` `--minwordmatches`,
//! `--weak-id` `--weak_id`, `--max-diffs` `--maxdiffs`, and so on.

use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// xorshift64, seeded as the generator of the expectations.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Rng {
        Rng(seed.wrapping_mul(2_654_435_761).wrapping_add(1))
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    fn base(&mut self) -> char {
        b"ACGT"[(self.next() >> 62) as usize] as char
    }

    fn seq(&mut self, n: usize) -> String {
        (0..n).map(|_| self.base()).collect()
    }

    /// `k` random substitutions, deletions, and insertions.
    fn edit(&mut self, seq: &str, k: usize) -> String {
        let mut x: Vec<char> = seq.chars().collect();
        for _ in 0..k {
            let op = self.below(3);
            let p = self.below(x.len() as u64) as usize;
            match op {
                0 => x[p] = self.base(),
                1 => {
                    x.remove(p);
                }
                _ => {
                    let c = self.base();
                    x.insert(p, c);
                }
            }
        }
        x.into_iter().collect()
    }
}

/// the reverse complement of `seq`, case kept.
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
            c => c,
        })
        .collect()
}

fn fasta(records: &[(&str, String)]) -> String {
    records
        .iter()
        .map(|(h, s)| format!(">{h}\n{s}\n"))
        .collect()
}

/// The inputs: a database of 32 sequences, most of them variants of one,
/// in mixed case, with Ns, one too short for 12 words, and the queries;
/// long sequences, and longer ones still, for the linear-memory aligner.
struct Inputs {
    db: String,
    q: String,
    ldb: String,
    lq: String,
    gdb: String,
    gq: String,
}

fn inputs() -> Inputs {
    let mut r = Rng::new(7);
    let base = r.seq(150);
    let mut db = vec![
        ("t0;size=4", base.clone()),
        ("t1", r.edit(&base, 2)),
        ("t2", r.edit(&base, 5)),
        ("t3", r.edit(&base, 10)),
        ("t4", r.edit(&base, 20)),
        ("t5", base[20..130].to_string()),
        ("t6", format!("{}{}", &base[..75], r.seq(75))),
        ("t7", r.seq(150)),
        ("t8;size=2", r.edit(&base, 3)),
        ("t9", r.edit(&base, 3).to_lowercase()),
        ("t10", base[..14].to_string()),
        ("t11", format!("{}NNNN{}", &base[..70], &base[74..])),
    ];
    let s = r.edit(&base, 1);
    db.push(("t12", format!("{}{}", s[..40].to_lowercase(), &s[40..])));
    let names: Vec<String> = (13..30).map(|i| format!("t{i}")).collect();
    for (i, name) in names.iter().enumerate() {
        db.push((name, r.edit(&base, 4 + i)));
    }
    db.push(("t30", r.seq(150)));
    db.push(("t31", r.seq(60)));
    let s = r.edit(&base, 6);
    let q = [
        ("q0", base.clone()),
        ("q1;size=3", r.edit(&base, 4)),
        ("q2", rc(&r.edit(&base, 2))),
        ("q3", r.seq(150)),
        ("q4", String::new()),
        ("q5", base[30..45].to_string()),
        ("q6", format!("{}{}", s[..75].to_lowercase(), &s[75..])),
        ("q7", format!("{}NN{}", &base[..60], &base[62..])),
        ("q8", format!("{}{}", &base[..75], r.seq(75))),
        ("q9", r.edit(&base, 15)),
        ("t0;size=4", r.edit(&base, 1)),
        ("q11", rc(&base)),
    ];
    let l0 = r.seq(3000);
    let ldb = [
        ("l0", l0.clone()),
        ("l1", r.edit(&l0, 30)),
        ("l2", r.seq(2990)),
    ];
    let lq = [("lq0", r.edit(&l0, 10)), ("lq1", rc(&r.edit(&l0, 5)))];
    let g0 = r.seq(5100);
    let gdb = [("g0", g0.clone()), ("g1", r.seq(5100))];
    let gq = [("gq0", r.edit(&g0, 20))];
    Inputs {
        db: fasta(&db),
        q: fasta(&q),
        ldb: fasta(&ldb),
        lq: fasta(&lq),
        gdb: fasta(&gdb),
        gq: fasta(&gq),
    }
}

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-dbsrch-{}-{name}", std::process::id()))
}

/// Runs dbsrch with `args` on the queries `q` and the database `db`.
fn run(name: &str, args: &[&str], q: &str, db: &str) -> Output {
    let (qp, dbp) = (
        temp_path(&format!("{name}-q")),
        temp_path(&format!("{name}-db")),
    );
    std::fs::write(&qp, q).unwrap();
    std::fs::write(&dbp, db).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dbsrch"))
        .arg("-i")
        .arg(&qp)
        .arg("--db")
        .arg(&dbp)
        .args(args)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let _ = std::fs::remove_file(&qp);
    let _ = std::fs::remove_file(&dbp);
    output
}

/// The thread counts every search is run with. Its hits must not change.
const THREADS: [&str; 4] = ["1", "2", "3", "8"];

/// Runs dbsrch on each of `THREADS`, which must succeed and write the same
/// hits each time, and returns its hits.
fn dbsrch(name: &str, args: &[&str], q: &str, db: &str) -> String {
    let mut first: Option<String> = None;
    for t in THREADS {
        let with_threads: Vec<&str> = args.iter().copied().chain(["-t", t]).collect();
        let output = run(name, &with_threads, q, db);
        assert!(
            output.status.success(),
            "dbsrch -t {t} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let hits = String::from_utf8(output.stdout).unwrap();
        match &first {
            None => first = Some(hits),
            Some(one) => assert_eq!(&hits, one, "-t {t} against -t 1, {args:?}"),
        }
    }
    first.unwrap()
}

/// Runs dbsrch with `args` at -t 1 and at each count of `counts`, and checks
/// that the output, the messages, and the exit code are the same.
fn same_on_threads(name: &str, args: &[&str], q: &str, db: &str, counts: &[&str]) -> Output {
    let with = |t: &str| -> Vec<String> {
        args.iter()
            .map(|a| a.to_string())
            .chain(["-t".to_string(), t.to_string()])
            .collect()
    };
    let one_args = with("1");
    let one = run(
        name,
        &one_args.iter().map(String::as_str).collect::<Vec<_>>(),
        q,
        db,
    );
    for &t in counts {
        let these = with(t);
        let other = run(
            name,
            &these.iter().map(String::as_str).collect::<Vec<_>>(),
            q,
            db,
        );
        assert_eq!(other.status.code(), one.status.code(), "-t {t}, {args:?}");
        assert!(
            other.stdout == one.stdout,
            "-t {t} writes other hits, {args:?}"
        );
        assert_eq!(
            String::from_utf8_lossy(&other.stderr),
            String::from_utf8_lossy(&one.stderr),
            "-t {t}, {args:?}"
        );
    }
    one
}

/// Many queries, more than a batch of every thread count: variants of the
/// database's sequences on either strand, unrelated ones, short and empty
/// ones, some with abundances, in mixed case.
fn many_queries(n: usize) -> String {
    let i = inputs();
    let targets: Vec<&str> = i.db.lines().skip(1).step_by(2).collect();
    let mut r = Rng::new(31);
    let mut records = Vec::with_capacity(n);
    for k in 0..n {
        let t = targets[r.below(targets.len() as u64) as usize];
        let seq = match k % 7 {
            0 => r.seq(150),
            1 => rc(&r.edit(t, 3)),
            2 => String::new(),
            3 => r.edit(t, 1).to_lowercase(),
            4 => t[..t.len().min(20)].to_string(),
            _ => r.edit(t, (k % 11) + 1),
        };
        let label = if k % 5 == 0 {
            format!("m{k};size={}", k % 9 + 1)
        } else {
            format!("m{k}")
        };
        records.push(format!(">{label}\n{seq}\n"));
    }
    records.concat()
}

#[test]
fn writes_the_same_hits_on_any_thread_count() {
    let i = inputs();
    let q = many_queries(1500);
    for args in [
        &["--id", "0.8"][..],
        &[
            "--id",
            "0.7",
            "--strand",
            "both",
            "--maxaccepts",
            "3",
            "--maxrejects",
            "4",
        ],
        &[
            "--id",
            "0.9",
            "--maxaccepts",
            "0",
            "--maxrejects",
            "0",
            "--top-hits-only",
        ],
        &[
            "--id",
            "0.6",
            "--weak-id",
            "0.5",
            "--maxhits",
            "2",
            "--strand",
            "both",
        ],
    ] {
        let one = same_on_threads("many", args, &q, &i.db, &["2", "3", "4", "7", "16", "0"]);
        assert!(one.status.success());
        // a hit or an N record for every query
        let queries: std::collections::BTreeSet<&str> = std::str::from_utf8(&one.stdout)
            .unwrap()
            .lines()
            .map(|l| l.split('\t').nth(8).unwrap())
            .collect();
        assert_eq!(queries.len(), 1500, "{args:?}");
    }
}

#[test]
fn stops_at_an_invalid_query_as_on_one_thread() {
    let i = inputs();
    let q = many_queries(1200);
    // an invalid abundance, and a character that cannot be read, after
    // enough hits to fill the output's buffer several times over, and as
    // the first query
    for bad in [">bad;size=0\nACGT\n", ">bad\nAC\u{1}GT\n"] {
        for input in [format!("{q}{bad}{q}"), format!("{bad}{q}")] {
            let one = same_on_threads("bad", &["--id", "0.8"], &input, &i.db, &["2", "3", "8"]);
            assert_eq!(one.status.code(), Some(1));
        }
    }
}

#[test]
fn takes_the_thread_count() {
    let i = inputs();
    let expected = dbsrch("t", &["--id", "0.9"], &i.q, &i.db);
    for t in [&["-t", "0"][..], &["--threads", "5"], &["-t", "1025"], &[]] {
        let args: Vec<&str> = ["--id", "0.9"]
            .iter()
            .copied()
            .chain(t.iter().copied())
            .collect();
        let output = run("t", &args, &i.q, &i.db);
        assert!(output.status.success(), "{t:?}");
        assert_eq!(String::from_utf8(output.stdout).unwrap(), expected, "{t:?}");
    }
    for t in ["-1", "x", "1.5", ""] {
        assert_eq!(
            run("t", &["--id", "0.9", "-t", t], &i.q, &i.db)
                .status
                .code(),
            Some(2),
            "{t:?}"
        );
    }
    assert_eq!(
        same_on_threads("empty", &["--id", "0.9"], "", &i.db, &["2", "8"]).stdout,
        b""
    );
}

#[test]
fn finds_the_best_hit() {
    let i = inputs();
    assert_eq!(
        dbsrch("id90", &["--id", "0.9"], &i.q, &i.db),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch("id97", &["--id", "0.97"], &i.q, &i.db),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
N	*	*	*	.	*	*	*	q6	*
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch("id0", &["--id", "0"], &i.q, &i.db),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	10	0	0.0	+	0	0	14I	q4	t10
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

/// the targets sharing the most words first, the shorter, then the first,
/// of ties
#[test]
fn accepts_up_to_maxaccepts() {
    let i = inputs();
    assert_eq!(
        dbsrch("acc5", &["--id", "0.8", "--maxaccepts", "5"], &i.q, &i.db),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	13	147	88.8	+	0	0	15MD16MI12M2I58MI13MD10MI15MD5M	q9	t13
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "acc3rej1",
            &["--id", "0.9", "--maxaccepts", "3", "--maxrejects", "1"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
N	*	*	*	.	*	*	*	q8	*
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "acc20rej2",
            &["--id", "0.95", "--maxaccepts", "20", "--maxrejects", "2"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q0	t5
H	11	150	100.0	+	0	0	=	q0	t11
H	12	150	99.3	+	0	0	68MI82M	q0	t12
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q0	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q0	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q0	t3
H	18	150	96.0	+	0	0	104MI45MD	q0	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q0	t2
H	19	150	95.4	+	0	0	D17MI20M2D21MD18MI70M	q0	t19
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	8	15	100.0	+	0	0	31I15M107I	q5	t8;size=2
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	16	15	100.0	+	0	0	29I15M107I	q5	t16
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	21	15	100.0	+	0	0	30I15M106I	q5	t21
H	29	15	100.0	+	0	0	31I15M109I	q5	t29
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q7	t5
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q7	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q7	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q7	t3
H	18	150	96.0	+	0	0	104MI45MD	q7	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q7	t2
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	16	151	96.7	+	0	0	6MD63MD32MI48MI	t0;size=4	t16
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "acc8rej0",
            &["--id", "0.95", "--maxaccepts", "8", "--maxrejects", "0"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	2	150	95.4	+	0	0	29M2I86MI35M	q0	t2
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	12	151	96.7	+	0	0	29MD9MD30MI55MI26M	q1;size=3	t12
H	5	151	96.4	+	0	0	20D9MD9MD85MI6M20D	q1;size=3	t5
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	8	151	95.5	+	0	0	21MI8MD9MD23MI40MI22MI26M	q1;size=3	t8;size=2
H	13	151	95.4	+	0	0	15MD13MD9MD32MI32MD20MI26M	q1;size=3	t13
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	14	151	96.7	+	0	0	5MI65MD41MD38M	t0;size=4	t14
H	16	151	96.7	+	0	0	6MD63MD32MI48MI	t0;size=4	t16
H	17	151	96.0	+	0	0	27MD42MD57MD22M	t0;size=4	t17
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

/// --maxaccepts 0 and --maxrejects 0 take every target
#[test]
fn takes_every_target_with_zero() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "all",
            &["--id", "0.95", "--maxaccepts", "0", "--maxrejects", "0"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q0	t5
H	10	150	100.0	+	0	0	14M136D	q0	t10
H	11	150	100.0	+	0	0	=	q0	t11
H	12	150	99.3	+	0	0	68MI82M	q0	t12
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q0	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q0	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q0	t3
H	18	150	96.0	+	0	0	104MI45MD	q0	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q0	t2
H	19	150	95.4	+	0	0	D17MI20M2D21MD18MI70M	q0	t19
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	12	151	96.7	+	0	0	29MD9MD30MI55MI26M	q1;size=3	t12
H	5	151	96.4	+	0	0	20D9MD9MD85MI6M20D	q1;size=3	t5
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	8	151	95.5	+	0	0	21MI8MD9MD23MI40MI22MI26M	q1;size=3	t8;size=2
H	13	151	95.4	+	0	0	15MD13MD9MD32MI32MD20MI26M	q1;size=3	t13
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	8	15	100.0	+	0	0	31I15M107I	q5	t8;size=2
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	16	15	100.0	+	0	0	29I15M107I	q5	t16
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	21	15	100.0	+	0	0	30I15M106I	q5	t21
H	29	15	100.0	+	0	0	31I15M109I	q5	t29
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q7	t5
H	10	150	100.0	+	0	0	14M136D	q7	t10
H	11	150	100.0	+	0	0	=	q7	t11
H	12	150	99.3	+	0	0	68MI82M	q7	t12
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q7	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q7	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q7	t3
H	18	150	96.0	+	0	0	104MI45MD	q7	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q7	t2
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	10	151	100.0	+	0	0	14M137D	t0;size=4	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	5	151	99.1	+	0	0	20D50MD60M20D	t0;size=4	t5
H	12	151	98.7	+	0	0	151M	t0;size=4	t12
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	14	151	96.7	+	0	0	5MI65MD41MD38M	t0;size=4	t14
H	16	151	96.7	+	0	0	6MD63MD32MI48MI	t0;size=4	t16
H	17	151	96.0	+	0	0	27MD42MD57MD22M	t0;size=4	t17
H	3	151	95.4	+	0	0	I6MD28M2D33MD65MD14M	t0;size=4	t3
H	18	151	95.4	+	0	0	68MD36MI45MD	t0;size=4	t18
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn searches_both_strands() {
    let i = inputs();
    assert_eq!(
        dbsrch("both", &["--id", "0.9", "--strand", "both"], &i.q, &i.db),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	0	152	98.7	-	0	0	13MD63MD74M	q2	t0;size=4
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	0	150	100.0	-	0	0	=	q11	t0;size=4
"##
    );
    assert_eq!(
        dbsrch(
            "bothall",
            &[
                "--id",
                "0.97",
                "--strand",
                "both",
                "--maxaccepts",
                "0",
                "--maxrejects",
                "0"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q0	t5
H	10	150	100.0	+	0	0	14M136D	q0	t10
H	11	150	100.0	+	0	0	=	q0	t11
H	12	150	99.3	+	0	0	68MI82M	q0	t12
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q0	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	5	152	99.1	-	0	0	21D56MD54M20D	q2	t5
H	0	152	98.7	-	0	0	13MD63MD74M	q2	t0;size=4
H	11	152	98.7	-	0	0	13MD63MD74M	q2	t11
H	1	152	97.4	-	0	0	13MD63MD68MD5M	q2	t1
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	8	15	100.0	+	0	0	31I15M107I	q5	t8;size=2
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	16	15	100.0	+	0	0	29I15M107I	q5	t16
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	21	15	100.0	+	0	0	30I15M106I	q5	t21
H	29	15	100.0	+	0	0	31I15M109I	q5	t29
N	*	*	*	.	*	*	*	q6	*
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q7	t5
H	10	150	100.0	+	0	0	14M136D	q7	t10
H	11	150	100.0	+	0	0	=	q7	t11
H	12	150	99.3	+	0	0	68MI82M	q7	t12
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q7	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	10	151	100.0	+	0	0	14M137D	t0;size=4	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	5	151	99.1	+	0	0	20D50MD60M20D	t0;size=4	t5
H	12	151	98.7	+	0	0	151M	t0;size=4	t12
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	0	150	100.0	-	0	0	=	q11	t0;size=4
H	5	150	100.0	-	0	0	20D110M20D	q11	t5
H	10	150	100.0	-	0	0	14M136D	q11	t10
H	11	150	100.0	-	0	0	=	q11	t11
H	12	150	99.3	-	0	0	68MI82M	q11	t12
H	1	150	98.7	-	0	0	144MD5M	q11	t1
H	8	150	98.0	-	0	0	21MI40MI40MI49M	q11	t8;size=2
H	13	150	98.0	-	0	0	15MD54MI32MD47M	q11	t13
H	15	150	97.4	-	0	0	I40MI30MI35MI45M	q11	t15
H	14	150	97.4	-	0	0	5MI106MD38M	q11	t14
H	16	150	97.4	-	0	0	6MD95MI48MI	q11	t16
"##
    );
}

/// weak hits count as rejected, and are written after the accepted ones
#[test]
fn writes_weak_hits() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "weak",
            &["--id", "0.97", "--weak-id", "0.8", "--maxaccepts", "4"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	12	151	96.7	+	0	0	29MD9MD30MI55MI26M	q1;size=3	t12
H	5	151	96.4	+	0	0	20D9MD9MD85MI6M20D	q1;size=3	t5
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	8	151	95.5	+	0	0	21MI8MD9MD23MI40MI22MI26M	q1;size=3	t8;size=2
H	13	151	95.4	+	0	0	15MD13MD9MD32MI32MD20MI26M	q1;size=3	t13
H	15	151	94.8	+	0	0	I29MD42MI35MI18MI26M	q1;size=3	t15
H	14	151	94.8	+	0	0	5MI24MD9MD73MD11MI26M	q1;size=3	t14
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	3	151	93.4	+	0	0	I6MD22MD7M3D85MI11MD14M	q1;size=3	t3
H	18	151	93.4	+	0	0	29MD9MD66MI19MI25MD	q1;size=3	t18
H	17	151	93.4	+	0	0	28M2D9MD111M	q1;size=3	t17
H	2	151	92.2	+	0	0	117MI8MI26M	q1;size=3	t2
H	23	151	90.3	+	0	0	12MI17MD9MI15MD17MD34MI44M	q1;size=3	t23
H	24	151	89.6	+	0	0	29MD10M2D16MI8MI18MI23MD43M	q1;size=3	t24
H	25	151	89.6	+	0	0	30MI9M3D26MI18MD25MI18MD20M	q1;size=3	t25
H	21	151	89.5	+	0	0	29MD9MD56MI24MI31M	q1;size=3	t21
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
H	22	151	88.8	+	0	0	6MD22MD5M2D24MD63MI17MD8M	q1;size=3	t22
H	27	151	87.2	+	0	0	17MD11M2D8MD11M2D27MI11MI10MI13MD10MI13MI13M	q1;size=3	t27
H	20	151	86.3	+	0	0	25M2D12MD82M2I29M	q1;size=3	t20
H	26	151	86.1	+	0	0	38MD47MD17MD20MI10MD14MD	q1;size=3	t26
H	4	151	84.9	+	0	0	15MD23MD11M2D35M2I13MD23MI24M2D	q1;size=3	t4
H	29	151	83.4	+	0	0	10MI16MD12MD82M2I4M2I12MI13M	q1;size=3	t29
H	28	151	80.3	+	0	0	28M3D8MD21MI9M2I21MI27MI7MI20MD5M	q1;size=3	t28
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	1	152	94.8	+	0	0	27MI14MD25MD56MD21MD5M	q6	t1
H	5	152	94.7	+	0	0	20D7MI14MD25MD56MD7M20D	q6	t5
H	14	152	93.5	+	0	0	5MI22MI14MD25MD44MD11MD27M	q6	t14
H	16	152	93.5	+	0	0	6MD20MI14MD25MD35MI21MD27MI	q6	t16
H	13	152	93.5	+	0	0	15MD11MI14MD61MD20MD27M	q6	t13
H	15	152	92.9	+	0	0	I27MI79MI18MD27M	q6	t15
H	18	152	92.8	+	0	0	26MI15MD25MD37MI19MD26MD	q6	t18
H	17	152	92.8	+	0	0	41MD25MD54M2D28M	q6	t17
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	3	152	90.2	+	0	0	I6MD20MI8M3D29MD56MD12MD14M	q6	t3
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	2	152	89.0	+	0	0	30M3I11MD25MD84M	q6	t2
H	24	152	88.4	+	0	0	27MI14M2D13MI27MI23MD15M2D28M	q6	t24
H	22	152	88.2	+	0	0	6MD34MD17MD7MD56MD18MD8M	q6	t22
H	25	152	87.2	+	0	0	27M3I13M3D42MD25MI11M2D7MD20M	q6	t25
H	21	152	87.0	+	0	0	27MI14MD25MD27MI23MD33M	q6	t21
H	20	152	86.8	+	0	0	41MD25MD84M	q6	t20
H	23	152	86.5	+	0	0	12MI15MI13MI12MD10M2D41MI16M2D28M	q6	t23
H	27	152	86.5	+	0	0	17MD27M3D13MD17MI11MI10MI13MD10MD14MI13M	q6	t27
H	26	152	85.6	+	0	0	25M2I11MD30MD17MD17MD20MD11MD14MD	q6	t26
H	4	152	82.5	+	0	0	15MD10M2I19M3D19MD19M2I13MD23MD25M2D	q6	t4
H	29	152	82.2	+	0	0	10MI15MI16MD27MD51M2I18MI13M	q6	t29
H	28	152	80.0	+	0	0	29MD11MD17MI14MI17MI56MD5M	q6	t28
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	12	147	88.9	+	0	0	32MI12M2I57M2I14MD10MI15MD5M	q9	t12
H	13	147	88.8	+	0	0	15MD16MI12M2I58MI13MD10MI15MD5M	q9	t13
H	5	147	88.4	+	0	0	20D12MI12M2I14MD42M2I14MD11M20D	q9	t5
H	24	147	87.5	+	0	0	32MI15MI6MI27MI20MI11MD29MD5M	q9	t24
H	17	147	87.4	+	0	0	44M2I14MD42M2I14MD25MD5M	q9	t17
H	19	147	86.8	+	0	0	D17MI14MI26M2D18MI23M2I14MD10MI15MD5M	q9	t19
H	16	147	86.4	+	0	0	6MD25MI12M2I14MD41M3I15MD10MI21M	q9	t16
H	18	147	86.2	+	0	0	32MI12M2I14MD42M3I14MD10MI19M2D	q9	t18
H	3	147	86.1	+	0	0	I6MD27MD9M2I14MD42M2I14MD25MD5M	q9	t3
H	2	147	85.8	+	0	0	29M3I15M2I14MD42M2I25MI15MD5M	q9	t2
H	22	147	85.3	+	0	0	6MD95M2I13MD10MI12M2D7M	q9	t22
H	15	147	85.2	+	0	0	I32MI12M3I57M3I14MD10MI15MD5M	q9	t15
H	14	147	84.2	+	0	0	5MI27MI12M2I14MD67MI15MD5M	q9	t14
H	20	147	83.6	+	0	0	44M2I14MD42M2I24MI16MD5M	q9	t20
H	21	147	80.4	+	0	0	32MI12M2I14MD33MI6M2I43MD5M	q9	t21
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "weakall",
            &[
                "--id",
                "0.99",
                "--weak-id",
                "0.9",
                "--maxaccepts",
                "0",
                "--maxrejects",
                "0"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q0	t5
H	10	150	100.0	+	0	0	14M136D	q0	t10
H	11	150	100.0	+	0	0	=	q0	t11
H	12	150	99.3	+	0	0	68MI82M	q0	t12
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q0	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q0	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q0	t3
H	18	150	96.0	+	0	0	104MI45MD	q0	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q0	t2
H	19	150	95.4	+	0	0	D17MI20M2D21MD18MI70M	q0	t19
H	22	150	93.3	+	0	0	6MD25MD26MD81MD8M	q0	t22
H	24	150	92.8	+	0	0	42MD13MI8MI18MI23MD19MD24M	q0	t24
H	20	150	92.1	+	0	0	25MD93MI31M	q0	t20
H	23	150	91.6	+	0	0	12MI26M2I14MD17MD34MI20MD24M	q0	t23
H	27	150	90.9	+	0	0	17MD13MD17M2D27MI11MI10MI13MD24MI13M	q0	t27
H	21	150	90.7	+	0	0	94MI56M	q0	t21
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	12	151	96.7	+	0	0	29MD9MD30MI55MI26M	q1;size=3	t12
H	5	151	96.4	+	0	0	20D9MD9MD85MI6M20D	q1;size=3	t5
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	8	151	95.5	+	0	0	21MI8MD9MD23MI40MI22MI26M	q1;size=3	t8;size=2
H	13	151	95.4	+	0	0	15MD13MD9MD32MI32MD20MI26M	q1;size=3	t13
H	15	151	94.8	+	0	0	I29MD42MI35MI18MI26M	q1;size=3	t15
H	14	151	94.8	+	0	0	5MI24MD9MD73MD11MI26M	q1;size=3	t14
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	3	151	93.4	+	0	0	I6MD22MD7M3D85MI11MD14M	q1;size=3	t3
H	18	151	93.4	+	0	0	29MD9MD66MI19MI25MD	q1;size=3	t18
H	17	151	93.4	+	0	0	28M2D9MD111M	q1;size=3	t17
H	2	151	92.2	+	0	0	117MI8MI26M	q1;size=3	t2
H	23	151	90.3	+	0	0	12MI17MD9MI15MD17MD34MI44M	q1;size=3	t23
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	8	15	100.0	+	0	0	31I15M107I	q5	t8;size=2
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	16	15	100.0	+	0	0	29I15M107I	q5	t16
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	21	15	100.0	+	0	0	30I15M106I	q5	t21
H	29	15	100.0	+	0	0	31I15M109I	q5	t29
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	1	152	94.8	+	0	0	27MI14MD25MD56MD21MD5M	q6	t1
H	5	152	94.7	+	0	0	20D7MI14MD25MD56MD7M20D	q6	t5
H	14	152	93.5	+	0	0	5MI22MI14MD25MD44MD11MD27M	q6	t14
H	16	152	93.5	+	0	0	6MD20MI14MD25MD35MI21MD27MI	q6	t16
H	13	152	93.5	+	0	0	15MD11MI14MD61MD20MD27M	q6	t13
H	15	152	92.9	+	0	0	I27MI79MI18MD27M	q6	t15
H	18	152	92.8	+	0	0	26MI15MD25MD37MI19MD26MD	q6	t18
H	17	152	92.8	+	0	0	41MD25MD54M2D28M	q6	t17
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	3	152	90.2	+	0	0	I6MD20MI8M3D29MD56MD12MD14M	q6	t3
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q7	t5
H	10	150	100.0	+	0	0	14M136D	q7	t10
H	11	150	100.0	+	0	0	=	q7	t11
H	12	150	99.3	+	0	0	68MI82M	q7	t12
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q7	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q7	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q7	t3
H	18	150	96.0	+	0	0	104MI45MD	q7	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q7	t2
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	22	150	93.3	+	0	0	6MD25MD26MD81MD8M	q7	t22
H	24	150	92.8	+	0	0	42MD13MI8MI18MI23MD19MD24M	q7	t24
H	20	150	92.1	+	0	0	25MD93MI31M	q7	t20
H	23	150	91.6	+	0	0	12MI26M2I14MD17MD34MI20MD24M	q7	t23
H	27	150	90.9	+	0	0	17MD13MD17M2D27MI11MI10MI13MD24MI13M	q7	t27
H	21	150	90.7	+	0	0	94MI56M	q7	t21
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	10	151	100.0	+	0	0	14M137D	t0;size=4	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	5	151	99.1	+	0	0	20D50MD60M20D	t0;size=4	t5
H	12	151	98.7	+	0	0	151M	t0;size=4	t12
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	14	151	96.7	+	0	0	5MI65MD41MD38M	t0;size=4	t14
H	16	151	96.7	+	0	0	6MD63MD32MI48MI	t0;size=4	t16
H	17	151	96.0	+	0	0	27MD42MD57MD22M	t0;size=4	t17
H	3	151	95.4	+	0	0	I6MD28M2D33MD65MD14M	t0;size=4	t3
H	18	151	95.4	+	0	0	68MD36MI45MD	t0;size=4	t18
H	2	151	94.8	+	0	0	29M2I41MD45MI35M	t0;size=4	t2
H	8	151	94.8	+	0	0	21MI81MI49M	t0;size=4	t8;size=2
H	19	151	94.7	+	0	0	D17MI20M2D21MD8MD10MI70M	t0;size=4	t19
H	20	151	93.4	+	0	0	25MD49MD44MI31M	t0;size=4	t20
H	22	151	92.7	+	0	0	6MD25MD26MD10MD71MD8M	t0;size=4	t22
H	27	151	92.3	+	0	0	17MD13MD17M2D10MD17MI11MI10MI13MD24MI13M	t0;size=4	t27
H	23	151	91.0	+	0	0	12MI26M2I14MD17M2D34MI20MD24M	t0;size=4	t23
H	24	151	90.8	+	0	0	42MD13MI27MI23MD19MD24M	t0;size=4	t24
H	21	151	90.1	+	0	0	70MD24MI56M	t0;size=4	t21
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn keeps_the_first_hits() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxhits",
            &["--id", "0.8", "--maxaccepts", "6", "--maxhits", "2"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "top",
            &["--id", "0.8", "--maxaccepts", "10", "--top-hits-only"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

/// shorter words make more candidates; few words a query or target lower
/// the words to share
#[test]
fn takes_the_word_length() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "w3",
            &["--id", "0.9", "--word-len", "3", "--maxaccepts", "3"],
            &i.q,
            &i.db
        ),
        r##"H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	27	150	90.9	+	0	0	17MD13MD17M2D27MI11MI10MI13MD24MI13M	q0	t27
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	17	151	93.4	+	0	0	28M2D9MD111M	q1;size=3	t17
H	2	151	92.2	+	0	0	117MI8MI26M	q1;size=3	t2
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	17	15	93.3	+	0	0	29I15M104I	q5	t17
H	1	152	94.8	+	0	0	27MI14MD25MD56MD21MD5M	q6	t1
H	17	152	92.8	+	0	0	41MD25MD54M2D28M	q6	t17
H	3	152	90.2	+	0	0	I6MD20MI8M3D29MD56MD12MD14M	q6	t3
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	27	150	90.9	+	0	0	17MD13MD17M2D27MI11MI10MI13MD24MI13M	q7	t27
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	27	151	92.3	+	0	0	17MD13MD17M2D10MD17MI11MI10MI13MD24MI13M	t0;size=4	t27
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "w5",
            &["--id", "0.9", "--word-len", "5", "--maxaccepts", "3"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	13	151	95.4	+	0	0	15MD13MD9MD32MI32MD20MI26M	q1;size=3	t13
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "w15",
            &["--id", "0.9", "--word-len", "15", "--maxaccepts", "3"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	2	151	92.2	+	0	0	117MI8MI26M	q1;size=3	t2
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
N	*	*	*	.	*	*	*	q8	*
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "mwm1",
            &[
                "--id",
                "0.95",
                "--min-word-matches",
                "1",
                "--maxaccepts",
                "0",
                "--maxrejects",
                "0"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q0	t5
H	10	150	100.0	+	0	0	14M136D	q0	t10
H	11	150	100.0	+	0	0	=	q0	t11
H	12	150	99.3	+	0	0	68MI82M	q0	t12
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q0	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q0	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q0	t3
H	18	150	96.0	+	0	0	104MI45MD	q0	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q0	t2
H	19	150	95.4	+	0	0	D17MI20M2D21MD18MI70M	q0	t19
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	12	151	96.7	+	0	0	29MD9MD30MI55MI26M	q1;size=3	t12
H	5	151	96.4	+	0	0	20D9MD9MD85MI6M20D	q1;size=3	t5
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	8	151	95.5	+	0	0	21MI8MD9MD23MI40MI22MI26M	q1;size=3	t8;size=2
H	13	151	95.4	+	0	0	15MD13MD9MD32MI32MD20MI26M	q1;size=3	t13
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	8	15	100.0	+	0	0	31I15M107I	q5	t8;size=2
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	16	15	100.0	+	0	0	29I15M107I	q5	t16
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	21	15	100.0	+	0	0	30I15M106I	q5	t21
H	29	15	100.0	+	0	0	31I15M109I	q5	t29
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q7	t5
H	10	150	100.0	+	0	0	14M136D	q7	t10
H	11	150	100.0	+	0	0	=	q7	t11
H	12	150	99.3	+	0	0	68MI82M	q7	t12
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q7	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q7	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q7	t3
H	18	150	96.0	+	0	0	104MI45MD	q7	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q7	t2
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	10	151	100.0	+	0	0	14M137D	t0;size=4	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	5	151	99.1	+	0	0	20D50MD60M20D	t0;size=4	t5
H	12	151	98.7	+	0	0	151M	t0;size=4	t12
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	14	151	96.7	+	0	0	5MI65MD41MD38M	t0;size=4	t14
H	16	151	96.7	+	0	0	6MD63MD32MI48MI	t0;size=4	t16
H	17	151	96.0	+	0	0	27MD42MD57MD22M	t0;size=4	t17
H	3	151	95.4	+	0	0	I6MD28M2D33MD65MD14M	t0;size=4	t3
H	18	151	95.4	+	0	0	68MD36MI45MD	t0;size=4	t18
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "mwm60",
            &[
                "--id",
                "0.95",
                "--min-word-matches",
                "60",
                "--maxaccepts",
                "0",
                "--maxrejects",
                "0"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q0	t5
H	10	150	100.0	+	0	0	14M136D	q0	t10
H	11	150	100.0	+	0	0	=	q0	t11
H	12	150	99.3	+	0	0	68MI82M	q0	t12
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q0	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q0	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q0	t3
H	18	150	96.0	+	0	0	104MI45MD	q0	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q0	t2
H	19	150	95.4	+	0	0	D17MI20M2D21MD18MI70M	q0	t19
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	12	151	96.7	+	0	0	29MD9MD30MI55MI26M	q1;size=3	t12
H	5	151	96.4	+	0	0	20D9MD9MD85MI6M20D	q1;size=3	t5
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	8	151	95.5	+	0	0	21MI8MD9MD23MI40MI22MI26M	q1;size=3	t8;size=2
H	13	151	95.4	+	0	0	15MD13MD9MD32MI32MD20MI26M	q1;size=3	t13
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	8	15	100.0	+	0	0	31I15M107I	q5	t8;size=2
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	16	15	100.0	+	0	0	29I15M107I	q5	t16
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	21	15	100.0	+	0	0	30I15M106I	q5	t21
H	29	15	100.0	+	0	0	31I15M109I	q5	t29
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q7	t5
H	10	150	100.0	+	0	0	14M136D	q7	t10
H	11	150	100.0	+	0	0	=	q7	t11
H	12	150	99.3	+	0	0	68MI82M	q7	t12
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q7	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	17	150	96.7	+	0	0	27MD99MD22M	q7	t17
H	3	150	96.0	+	0	0	I6MD28M2D98MD14M	q7	t3
H	18	150	96.0	+	0	0	104MI45MD	q7	t18
H	2	150	95.4	+	0	0	29M2I86MI35M	q7	t2
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	10	151	100.0	+	0	0	14M137D	t0;size=4	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	5	151	99.1	+	0	0	20D50MD60M20D	t0;size=4	t5
H	12	151	98.7	+	0	0	151M	t0;size=4	t12
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	14	151	96.7	+	0	0	5MI65MD41MD38M	t0;size=4	t14
H	16	151	96.7	+	0	0	6MD63MD32MI48MI	t0;size=4	t16
H	17	151	96.0	+	0	0	27MD42MD57MD22M	t0;size=4	t17
H	3	151	95.4	+	0	0	I6MD28M2D33MD65MD14M	t0;size=4	t3
H	18	151	95.4	+	0	0	68MD36MI45MD	t0;size=4	t18
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "mwm0",
            &[
                "--id",
                "0.97",
                "--min-word-matches",
                "0",
                "--maxaccepts",
                "0",
                "--maxrejects",
                "0"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q0	t5
H	10	150	100.0	+	0	0	14M136D	q0	t10
H	11	150	100.0	+	0	0	=	q0	t11
H	12	150	99.3	+	0	0	68MI82M	q0	t12
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q0	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q0	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	9	150	97.3	+	0	0	I150M	q0	t9
H	10	151	100.0	+	0	0	14M137D	q1;size=3	t10
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	6	15	100.0	+	0	0	30I15M105I	q5	t6
H	8	15	100.0	+	0	0	31I15M107I	q5	t8;size=2
H	9	15	100.0	+	0	0	31I15M105I	q5	t9
H	11	15	100.0	+	0	0	30I15M105I	q5	t11
H	12	15	100.0	+	0	0	30I15M106I	q5	t12
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	14	15	100.0	+	0	0	31I15M104I	q5	t14
H	16	15	100.0	+	0	0	29I15M107I	q5	t16
H	18	15	100.0	+	0	0	30I15M105I	q5	t18
H	20	15	100.0	+	0	0	29I15M106I	q5	t20
H	21	15	100.0	+	0	0	30I15M106I	q5	t21
H	29	15	100.0	+	0	0	31I15M109I	q5	t29
H	10	152	100.0	+	0	0	14M138D	q6	t10
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	5	150	100.0	+	0	0	20D110M20D	q7	t5
H	10	150	100.0	+	0	0	14M136D	q7	t10
H	11	150	100.0	+	0	0	=	q7	t11
H	12	150	99.3	+	0	0	68MI82M	q7	t12
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	15	150	97.4	+	0	0	I40MI30MI35MI45M	q7	t15
H	14	150	97.4	+	0	0	5MI106MD38M	q7	t14
H	16	150	97.4	+	0	0	6MD95MI48MI	q7	t16
H	9	150	97.3	+	0	0	I150M	q7	t9
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	10	151	100.0	+	0	0	14M137D	t0;size=4	t10
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	5	151	99.1	+	0	0	20D50MD60M20D	t0;size=4	t5
H	12	151	98.7	+	0	0	151M	t0;size=4	t12
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn takes_identity_definition_0() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "iddef0",
            &["--id", "0.9", "--iddef", "0", "--maxaccepts", "5"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	8	150	100.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	99.3	+	0	0	144MD5M	q0	t1
H	13	150	99.3	+	0	0	15MD54MI32MD47M	q0	t13
H	0	151	98.7	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	98.7	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	98.0	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	96.0	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	19	151	91.9	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	100.0	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	98.0	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	98.0	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	96.7	+	0	0	27MI14MD82MD27M	q6	t12
H	8	152	94.1	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	19	152	93.2	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	8	150	100.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	99.3	+	0	0	144MD5M	q7	t1
H	19	150	97.3	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	11	147	95.2	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	94.6	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	8	147	93.9	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	1	147	93.2	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	13	147	91.8	+	0	0	15MD16MI12M2I58MI13MD10MI15MD5M	q9	t13
H	0	151	100.0	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	100.0	+	0	0	70MD80M	t0;size=4	t11
H	15	151	99.3	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	99.3	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	99.3	+	0	0	15MD87MD47M	t0;size=4	t13
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn takes_identity_definition_1() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "iddef1",
            &["--id", "0.9", "--iddef", "1", "--maxaccepts", "5"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.2	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	2	151	92.2	+	0	0	117MI8MI26M	q1;size=3	t2
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	13	152	93.5	+	0	0	15MD11MI14MD61MD20MD27M	q6	t13
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	19	150	94.7	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
N	*	*	*	.	*	*	*	q8	*
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	15	151	97.4	+	0	0	I40MI66MI45M	t0;size=4	t15
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn takes_identity_definition_3() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "iddef3",
            &["--id", "0.9", "--iddef", "3", "--maxaccepts", "5"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.0	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.0	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	19	151	90.1	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	0	15	98.7	+	0	0	30I15M105I	q5	t0;size=4
H	1	15	98.7	+	0	0	30I15M104I	q5	t1
H	13	15	98.7	+	0	0	29I15M105I	q5	t13
H	4	15	98.6	+	0	0	30I15M102I	q5	t4
H	5	15	98.2	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	8	152	92.2	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	19	152	90.8	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	19	150	95.3	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	10	150	99.3	+	0	0	14M136D	q8	t10
H	11	147	92.7	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	92.0	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	91.3	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	5	147	91.2	+	0	0	20D12MI12M2I14MD42M2I14MD11M20D	q9	t5
H	8	147	90.8	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	15	151	97.4	+	0	0	I40MI66MI45M	t0;size=4	t15
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn takes_identity_definition_4() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "iddef4",
            &["--id", "0.9", "--iddef", "4", "--maxaccepts", "5"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.2	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	2	151	92.2	+	0	0	117MI8MI26M	q1;size=3	t2
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	13	152	93.5	+	0	0	15MD11MI14MD61MD20MD27M	q6	t13
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	19	150	94.7	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
N	*	*	*	.	*	*	*	q8	*
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
H	15	151	97.4	+	0	0	I40MI66MI45M	t0;size=4	t15
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxdiffs() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxdiffs",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-diffs",
                "2"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
N	*	*	*	.	*	*	*	q6	*
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxsubs() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxsubs",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-subs",
                "0"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
N	*	*	*	.	*	*	*	q6	*
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxgaps() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxgaps",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-gaps",
                "1"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
N	*	*	*	.	*	*	*	q6	*
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
N	*	*	*	.	*	*	*	q9	*
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_mincols() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "mincols",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--min-cols",
                "140"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_leftjust() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "leftjust",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--left-just"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_rightjust() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "rightjust",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--right-just"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_qcov() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "qcov",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--query-cov",
                "0.9"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_tcov() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "tcov",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--target-cov",
                "0.9"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxid() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxid",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-id",
                "0.97"
            ],
            &i.q,
            &i.db
        ),
        r##"N	*	*	*	.	*	*	*	q0	*
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	2	151	92.2	+	0	0	117MI8MI26M	q1;size=3	t2
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
N	*	*	*	.	*	*	*	q7	*
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
N	*	*	*	.	*	*	*	t0;size=4	*
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_mid() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "mid",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--mid",
                "99"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
N	*	*	*	.	*	*	*	q6	*
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_idprefix() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "idprefix",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--id-prefix",
                "5"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_idsuffix() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "idsuffix",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--id-suffix",
                "3"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
N	*	*	*	.	*	*	*	q8	*
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxqt() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxqt",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-qt",
                "0.95"
            ],
            &i.q,
            &i.db
        ),
        r##"N	*	*	*	.	*	*	*	q0	*
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
N	*	*	*	.	*	*	*	q6	*
N	*	*	*	.	*	*	*	q7	*
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
N	*	*	*	.	*	*	*	t0;size=4	*
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_minqt() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "minqt",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--min-qt",
                "1.05"
            ],
            &i.q,
            &i.db
        ),
        r##"N	*	*	*	.	*	*	*	q0	*
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
N	*	*	*	.	*	*	*	q6	*
N	*	*	*	.	*	*	*	q7	*
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
N	*	*	*	.	*	*	*	t0;size=4	*
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_minsl() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "minsl",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--min-sl",
                "0.9"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxsl() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxsl",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-sl",
                "0.95"
            ],
            &i.q,
            &i.db
        ),
        r##"N	*	*	*	.	*	*	*	q0	*
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
N	*	*	*	.	*	*	*	q6	*
N	*	*	*	.	*	*	*	q7	*
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
N	*	*	*	.	*	*	*	t0;size=4	*
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxqsize() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxqsize",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-qsize",
                "2"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
N	*	*	*	.	*	*	*	t0;size=4	*
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_mintsize() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "mintsize",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--min-tsize",
                "2"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_minsizer() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "minsizer",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--min-size-ratio",
                "1"
            ],
            &i.q,
            &i.db
        ),
        r##"H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	16	150	97.4	+	0	0	6MD95MI48MI	q0	t16
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	13	152	93.5	+	0	0	15MD11MI14MD61MD20MD27M	q6	t13
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q7	t13
H	19	150	95.4	+	0	0	D17MI20M2D20MD19MI70M	q7	t19
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	16	150	68.4	+	0	0	6MD86M3I56MD	q8	t16
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	13	147	88.8	+	0	0	15MD16MI12M2I58MI13MD10MI15MD5M	q9	t13
H	15	147	85.2	+	0	0	I32MI12M3I57M3I14MD10MI15MD5M	q9	t15
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_maxsizer() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "maxsizer",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--max-size-ratio",
                "0.5"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
N	*	*	*	.	*	*	*	t0;size=4	*
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_self() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "self",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--self"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	13	151	98.0	+	0	0	15MD87MD47M	t0;size=4	t13
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn applies_selfid() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "selfid",
            &[
                "--id",
                "0",
                "--maxaccepts",
                "4",
                "--maxrejects",
                "4",
                "--self-id"
            ],
            &i.q,
            &i.db
        ),
        r##"H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	13	150	98.0	+	0	0	15MD54MI32MD47M	q0	t13
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	89.6	+	0	0	D17MI9MI13M3D19M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	6	150	75.3	+	0	0	75M2I12M2I20MI14M4I20M9D	q8	t6
H	1	150	68.4	+	0	0	93M2I18M3D36M	q8	t1
H	0	150	68.0	+	0	0	93M2I55M2D	q8	t0;size=4
H	11	150	68.0	+	0	0	93M2I55M2D	q8	t11
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	1	147	89.5	+	0	0	32MI12M2I14MD42M2I14MD10MI14M2D5M	q9	t1
H	8	147	89.0	+	0	0	21MI11MI12M2I56M3I15MD10MI15MD5M	q9	t8;size=2
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn takes_the_scores() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "score",
            &[
                "--id",
                "0.8",
                "--maxaccepts",
                "4",
                "--match",
                "1",
                "--mismatch",
                "-2"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	1	150	97.3	+	0	0	149MD	q0	t1
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	94.7	+	0	0	29MD9MD85MI25MD	q1;size=3	t1
H	19	151	89.5	+	0	0	D38M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	90.1	+	0	0	67MD56MD27M	q6	t0;size=4
H	11	152	90.1	+	0	0	67MD56MD27M	q6	t11
H	12	152	89.5	+	0	0	124MD27M	q6	t12
H	19	152	83.6	+	0	0	D21MD39M2D17MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	1	150	97.3	+	0	0	149MD	q7	t1
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	11	147	84.0	+	0	0	32MI26MI43M2I45MD	q9	t11
H	8	147	83.7	+	0	0	21MI11MI12M2I56M3I46MD	q9	t8;size=2
H	0	147	83.3	+	0	0	32MI26MI43M2I45MD	q9	t0;size=4
H	1	147	81.2	+	0	0	32MI26MI43M2I44M2D	q9	t1
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	96.7	+	0	0	70MD79MD	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "gaps",
            &[
                "--id",
                "0.8",
                "--maxaccepts",
                "4",
                "--gapopen",
                "10I/1E",
                "--gapext",
                "1I/1E"
            ],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	19	151	92.8	+	0	0	D17MI11MD9M3D21MD18MI43MI26M	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	19	152	91.6	+	0	0	D17MI9MI13M3D18MD5MD13MI43MD27M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	98.7	+	0	0	144MD5M	q7	t1
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	11	147	92.8	+	0	0	32MI12MI3MI11MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	92.2	+	0	0	32MI12MI3MI11MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	8	147	91.6	+	0	0	21MI11MI12MI3MI53M2I2MI13MD10MI15MD5M	q9	t8;size=2
H	1	147	90.8	+	0	0	32MI12MI3MI11MD42M2I14MD10MI14M2D5M	q9	t1
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "nmis",
            &["--id", "0.9", "--maxaccepts", "4", "--n-mismatch"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	8	150	98.0	+	0	0	21MI40MI40MI49M	q0	t8;size=2
H	11	150	97.3	+	0	0	150M	q0	t11
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	1	151	96.1	+	0	0	29MD9MD85MI20MD5M	q1;size=3	t1
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
H	11	151	94.7	+	0	0	29MD9MD85MI26M	q1;size=3	t11
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	11	152	93.5	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	8	152	92.3	+	0	0	21MI6MI14MD60MI22MD27M	q6	t8;size=2
H	0	150	98.7	+	0	0	150M	q7	t0;size=4
H	1	150	97.3	+	0	0	144MD5M	q7	t1
H	8	150	96.7	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	11	150	96.0	+	0	0	150M	q7	t11
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
H	1	151	98.0	+	0	0	70MD74MD5M	t0;size=4	t1
H	11	151	96.7	+	0	0	70MD80M	t0;size=4	t11
N	*	*	*	.	*	*	*	q11	*
"##
    );
    assert_eq!(
        dbsrch(
            "inf",
            &["--id", "0", "--maxaccepts", "4", "--gapopen", "*I/2E"],
            &i.q,
            &i.db
        ),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	97.3	+	0	0	149MD	q0	t1
H	8	150	56.7	+	0	0	3I150M	q0	t8;size=2
H	1	151	71.8	+	0	0	2D149M	q1;size=3	t1
H	0	151	70.5	+	0	0	2D149MI	q1;size=3	t0;size=4
H	11	151	70.5	+	0	0	2D149MI	q1;size=3	t11
H	19	151	59.9	+	0	0	4D147MI	q1;size=3	t19
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
H	5	0	0.0	+	0	0	110I	q4	t5
H	10	0	0.0	+	0	0	14I	q4	t10
H	22	0	0.0	+	0	0	146I	q4	t22
H	31	0	0.0	+	0	0	60I	q4	t31
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	13	15	100.0	+	0	0	29I15M105I	q5	t13
H	12	152	78.8	+	0	0	151MD	q6	t12
H	0	152	53.3	+	0	0	D150MD	q6	t0;size=4
H	11	152	53.3	+	0	0	D150MD	q6	t11
H	19	152	50.7	+	0	0	4D148M	q6	t19
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	1	150	97.3	+	0	0	149MD	q7	t1
H	8	150	57.3	+	0	0	3I150M	q7	t8;size=2
H	0	150	62.7	+	0	0	150M	q8	t0;size=4
H	11	150	62.7	+	0	0	150M	q8	t11
H	1	150	62.4	+	0	0	149MD	q8	t1
H	6	150	60.7	+	0	0	150M	q8	t6
H	8	147	49.0	+	0	0	4I147M2I	q9	t8;size=2
H	11	147	46.9	+	0	0	3I147M	q9	t11
H	1	147	46.3	+	0	0	2I147M	q9	t1
H	0	147	44.9	+	0	0	3I147M	q9	t0;size=4
H	0	151	66.0	+	0	0	D150M	t0;size=4	t0;size=4
H	11	151	66.0	+	0	0	D150M	t0;size=4	t11
H	1	151	63.8	+	0	0	149M2D	t0;size=4	t1
H	15	151	60.3	+	0	0	2I151MI	t0;size=4	t15
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

/// scores beyond 16 bits fall back to the linear-memory aligner
#[test]
fn aligns_long_sequences() {
    let i = inputs();
    assert_eq!(
        dbsrch(
            "long",
            &["--id", "0.5", "--strand", "both", "--maxaccepts", "0"],
            &i.lq,
            &i.ldb
        ),
        r##"H	0	3005	99.7	+	0	0	29MD127MD926MI276MD15MD349MD31MD1246M	lq0	l0
H	1	3005	98.8	+	0	0	29MD49MD77MD565MD82MD125MI107MI19MD25MI131MD149MD10MD30MI264MD54MD24MD6MD32MD129MI281MI68MI15MI311MI85MD98MD155MD68M	lq0	l1
H	0	3002	99.8	-	0	0	1298MI534MD606MD359MD202M	lq1	l0
H	1	3002	98.9	-	0	0	78MD642MD82MD125MI107MI19MD157MD83MI105MI264MD78MD38MD46MD83MI281MI68MI15MI159MD152MI85MD98MD22MD133MD68M	lq1	l1
"##
    );
    assert_eq!(
        dbsrch(
            "long12",
            &[
                "--id",
                "0.5",
                "--strand",
                "both",
                "--maxaccepts",
                "0",
                "--match",
                "12"
            ],
            &i.lq,
            &i.ldb
        ),
        r##"H	0	3005	99.7	+	0	0	29MD127MD926MI276MD15MD349MD31MD1246M	lq0	l0
H	1	3005	98.8	+	0	0	29MD49MD77MD565MD82MD125MI107MI19MD25MI131MD149MD10MD30MI264MD54MD24MD6MD32MD130MI280MI68MI15MI314MI82MD98MD155MD68M	lq0	l1
H	0	3002	99.8	-	0	0	1300MI532MD606MD359MD202M	lq1	l0
H	1	3002	98.9	-	0	0	78MD642MD82MD125MI107MI19MD157MD85MI103MI264MD78MD38MD46MD84MI280MI68MI15MI159MD155MI82MD98MD22MD133MD68M	lq1	l1
"##
    );
}

#[test]
fn aligns_sequences_too_long_for_16_bits() {
    let i = inputs();
    assert_eq!(
        dbsrch("big", &["--id", "0.5", "--maxaccepts", "0"], &i.gq, &i.gdb),
        r##"H	0	5102	99.7	+	0	0	99MI119MD30MD21MI176MI556MI25MD323MD1139MD715MD81MD742MD753MI141MI174M	gq0	g0
"##
    );
}

#[test]
fn reads_udb_databases() {
    // a UDB file's word length wins, with a warning when it differs
    let i = inputs();
    let (db, udb) = (temp_path("udb-db"), temp_path("udb-udb"));
    std::fs::write(&db, &i.db).unwrap();
    let made = Command::new(env!("CARGO_BIN_EXE_mkudb"))
        .args([
            "-i",
            db.to_str().unwrap(),
            "-o",
            udb.to_str().unwrap(),
            "--word-len",
            "6",
        ])
        .output()
        .unwrap();
    assert!(made.status.success());
    let search = |args: &[&str]| {
        let q = temp_path(&format!("udb-q{}", args.len()));
        std::fs::write(&q, &i.q).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dbsrch"))
            .args(["-i", q.to_str().unwrap(), "--db", udb.to_str().unwrap()])
            .args(["--id", "0.9", "--maxaccepts", "3"])
            .args(args)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        let _ = std::fs::remove_file(&q);
        output
    };
    let warned = search(&[]);
    let quiet = search(&["--word-len", "6"]);
    let _ = std::fs::remove_file(&db);
    let _ = std::fs::remove_file(&udb);
    let expected = r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	11	150	100.0	+	0	0	=	q0	t11
H	1	150	98.7	+	0	0	144MD5M	q0	t1
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
H	11	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t11
H	16	151	94.8	+	0	0	6MD22MD9MD64MI21MI26MI	q1;size=3	t16
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	1	15	100.0	+	0	0	30I15M104I	q5	t1
H	4	15	100.0	+	0	0	30I15M102I	q5	t4
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	11	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t11
H	12	152	95.4	+	0	0	27MI14MD82MD27M	q6	t12
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	11	150	100.0	+	0	0	=	q7	t11
H	8	150	98.0	+	0	0	21MI39MI41MI49M	q7	t8;size=2
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	10	147	100.0	+	0	0	14M133D	q9	t10
H	11	147	91.5	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t11
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
H	11	151	99.3	+	0	0	70MD80M	t0;size=4	t11
H	15	151	98.0	+	0	0	I40MI66MI45M	t0;size=4	t15
N	*	*	*	.	*	*	*	q11	*
"##;
    assert!(warned.status.success());
    assert_eq!(String::from_utf8(warned.stdout).unwrap(), expected);
    assert!(String::from_utf8_lossy(&warned.stderr).contains("word length"));
    assert!(quiet.status.success());
    assert_eq!(String::from_utf8(quiet.stdout).unwrap(), expected);
    assert!(
        quiet.stderr.is_empty(),
        "{}",
        String::from_utf8_lossy(&quiet.stderr)
    );
    assert_eq!(
        dbsrch(
            "udb-fa",
            &["--id", "0.9", "--maxaccepts", "3", "--word-len", "6"],
            &i.q,
            &i.db
        ),
        expected
    );
}

/// FASTA records as FASTQ ones.
fn fastq(text: &str) -> String {
    text.split('>')
        .skip(1)
        .map(|record| {
            let (header, seq) = record.split_once('\n').unwrap();
            let seq = seq.trim_end();
            format!("@{header}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
        })
        .collect()
}

#[test]
fn reads_fastq() {
    let i = inputs();
    let args = ["--id", "0.9"];
    assert_eq!(
        dbsrch("fastq", &args, &fastq(&i.q), &fastq(&i.db)),
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn searches_empty_files() {
    let i = inputs();
    assert_eq!(dbsrch("empty-q", &["--id", "0.9"], "", &i.db), "");
    // without a database, no query has a hit
    assert_eq!(
        dbsrch("empty-db", &["--id", "0.9"], &i.q, ""),
        r##"N	*	*	*	.	*	*	*	q0	*
N	*	*	*	.	*	*	*	q1;size=3	*
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
N	*	*	*	.	*	*	*	q5	*
N	*	*	*	.	*	*	*	q6	*
N	*	*	*	.	*	*	*	q7	*
N	*	*	*	.	*	*	*	q8	*
N	*	*	*	.	*	*	*	q9	*
N	*	*	*	.	*	*	*	t0;size=4	*
N	*	*	*	.	*	*	*	q11	*
"##
    );
}

#[test]
fn checks_abundances() {
    let i = inputs();
    // the database is read before the queries
    let output = run("bad-db", &["--id", "0.9"], &i.q, ">a;size=0\nACGT\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    assert!(String::from_utf8_lossy(&output.stderr).contains("abundance"));
    let output = run(
        "bad-q",
        &["--id", "0.9"],
        &format!("{}>b;size=0\nACGT\n", i.q),
        &i.db,
    );
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("abundance"));
}

#[test]
fn rejects_invalid_arguments() {
    let i = inputs();
    for args in [
        &[][..],
        &["--id", "1.5"],
        &["--id", "0.9", "--strand", "minus"],
        &["--id", "0.9", "--word-len", "2"],
        &["--id", "0.9", "--word-len", "16"],
        &["--id", "0.9", "--maxhits", "0"],
        &["--id", "0.9", "--maxaccepts", "-1"],
        &["--id", "0.9", "--iddef", "5"],
    ] {
        assert_eq!(
            run("args", args, &i.q, &i.db).status.code(),
            Some(2),
            "{args:?}"
        );
    }
    let output = Command::new(env!("CARGO_BIN_EXE_dbsrch"))
        .args(["--id", "0.9"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    let missing = temp_path("missing");
    let output = Command::new(env!("CARGO_BIN_EXE_dbsrch"))
        .args(["--id", "0.9", "--db", missing.to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
}

#[test]
fn reads_stdin_and_writes_files() {
    let i = inputs();
    let (db, out) = (temp_path("io-db"), temp_path("io-out"));
    std::fs::write(&db, &i.db).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_dbsrch"))
        .args([
            "--id",
            "0.9",
            "--db",
            db.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    use std::io::Write;
    child
        .stdin
        .take()
        .unwrap()
        .write_all(i.q.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    let _ = std::fs::remove_file(&db);
    let _ = std::fs::remove_file(&out);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(
        written,
        r##"H	0	150	100.0	+	0	0	=	q0	t0;size=4
H	0	151	97.4	+	0	0	29MD9MD85MI26M	q1;size=3	t0;size=4
N	*	*	*	.	*	*	*	q2	*
N	*	*	*	.	*	*	*	q3	*
N	*	*	*	.	*	*	*	q4	*
H	5	15	100.0	+	0	0	10I15M85I	q5	t5
H	0	152	96.1	+	0	0	27MI14MD25MD56MD27M	q6	t0;size=4
H	0	150	100.0	+	0	0	=	q7	t0;size=4
H	10	150	100.0	+	0	0	14M136D	q8	t10
H	0	147	90.8	+	0	0	32MI12M2I14MD42M2I14MD10MI15MD5M	q9	t0;size=4
H	0	151	99.3	+	0	0	70MD80M	t0;size=4	t0;size=4
N	*	*	*	.	*	*	*	q11	*
"##
    );
}
