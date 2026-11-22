//! End-to-end tests of the allpr binary.
//!
//! Expected hits are those of `vsearch --allpairs_global --uc
//! --notrunclabels --qmask none --minseqlength 1` on the same input. Every
//! sequence is aligned with each one after it, globally, by vsearch's
//! 16-bit aligner (or its linear-memory one when a score overflows or the
//! sequences are too long). The hits accepted by vsearch's criteria come by
//! decreasing identity, then target. The options are vsearch's, in words:
//! `--accept-all` is `--acceptall`, `--top-hits-only` `--top_hits_only`,
//! `--weak-id` `--weak_id`, `--max-diffs` `--maxdiffs`, and so on.

use std::io::Write;
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

/// The inputs: ten related or random sequences, two of one label then a
/// third, three long ones, two too long for the 16-bit aligner, and empty
/// ones.
struct Inputs {
    main: String,
    selves: String,
    long: String,
    big: String,
}

fn inputs() -> Inputs {
    let mut r = Rng::new(5);
    let base = r.seq(80);
    let seqs = [
        ("s0;size=5", base.clone()),
        ("s1", r.edit(&base, 3)),
        ("s2", r.edit(&base, 8)),
        ("s3;size=2", base[10..70].to_string()),
        ("s4", format!("{}{}", r.edit(&base[..50], 2), r.seq(30))),
        ("s5", r.seq(80)),
        ("s6", format!("{}NNRY{}", &base[..40], &base[44..])),
        ("s7", r.edit(&base, 1).to_lowercase()),
        ("s8", r.edit(&base, 2)),
        ("s9", format!("{}{}", &base[5..], r.seq(12))),
    ];
    let main = seqs.iter().map(|(h, s)| format!(">{h}\n{s}\n")).collect();
    let selves = format!(
        ">a\n{base}\n>a\n{}\n>b\n{}\n",
        r.edit(&base, 2),
        r.edit(&base, 1)
    );
    let long_seq = r.seq(4000);
    let long = format!(
        ">l0\n{long_seq}\n>l1\n{}\n>l2\n{}\n",
        r.edit(&long_seq, 30),
        r.seq(3990)
    );
    let big_seq = r.seq(5100);
    let big = format!(">g0\n{big_seq}\n>g1\n{}\n", r.edit(&big_seq, 20));
    Inputs {
        main,
        selves,
        long,
        big,
    }
}

/// Runs allpr with `args` on `input`, fed on stdin, at the SIMD level
/// `level` or the processor's best.
fn run_at(level: Option<rsearch::simd::Level>, args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_allpr"))
        .envs(level.map(|l| (rsearch::simd::VARIABLE, l.name())))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("allpr should start");
    let mut stdin = child.stdin.take().unwrap();
    let input = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&input);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// The thread counts every run is made with. Its hits must not change.
const THREADS: [&str; 4] = ["1", "2", "3", "8"];

/// Runs allpr on each of `THREADS`, which must succeed and write the same
/// hits each time, and returns its hits.
fn allpr(args: &[&str], input: &str) -> String {
    let mut first: Option<String> = None;
    for t in THREADS {
        let with_threads: Vec<&str> = args.iter().copied().chain(["-t", t]).collect();
        let output = run(&with_threads, input);
        assert!(
            output.status.success(),
            "allpr -t {t} failed: {}",
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

/// Runs allpr with `args` at -t 1 and at each count of `counts`, and checks
/// that the hits, the messages, and the exit code are the same.
fn same_on_threads(args: &[&str], input: &str, counts: &[&str]) -> Output {
    let with = |t: &str| -> Vec<String> {
        args.iter()
            .map(|a| a.to_string())
            .chain(["-t".to_string(), t.to_string()])
            .collect()
    };
    let one_args = with("1");
    let one = run(
        &one_args.iter().map(String::as_str).collect::<Vec<_>>(),
        input,
    );
    for &t in counts {
        let these = with(t);
        let other = run(&these.iter().map(String::as_str).collect::<Vec<_>>(), input);
        assert_eq!(
            other.status.code(),
            one.status.code(),
            "-t {t}, {args:?}: {}",
            String::from_utf8_lossy(&other.stderr)
        );
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

/// Many sequences, more than a batch of every thread count: variants of a few
/// templates, unrelated, short, and empty ones, some with abundances.
fn many_sequences(n: usize) -> String {
    let mut r = Rng::new(23);
    let templates: Vec<String> = (0..4).map(|_| r.seq(90)).collect();
    let mut out = String::new();
    for k in 0..n {
        let t = &templates[k % templates.len()];
        let seq = match k % 8 {
            0 => r.seq(90),
            1 => String::new(),
            2 => t[..25].to_string(),
            3 => r.edit(t, 1).to_lowercase(),
            _ => r.edit(t, k % 9 + 1),
        };
        let label = if k % 5 == 0 {
            format!("m{k};size={}", k % 7 + 1)
        } else {
            format!("m{k}")
        };
        out.push_str(&format!(">{label}\n{seq}\n"));
    }
    out
}

#[test]
fn aligns_every_pair() {
    assert_eq!(
        allpr(&["--accept-all"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn takes_the_identity() {
    let main = inputs().main;
    assert_eq!(
        allpr(&["--id", "0.9"], &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
    assert_eq!(
        allpr(&["--id", "0.99"], &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
N	*	*	*	.	*	*	*	s1	*
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
N	*	*	*	.	*	*	*	s6	*
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
    // every hit passes --id 0
    assert_eq!(
        allpr(&["--id", "0"], &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn takes_identity_definition_0() {
    let args = ["--id", "0.9", "--iddef", "0"];
    assert_eq!(
        allpr(&args, &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	98.7	+	0	0	41MD38M	s0;size=5	s8
H	1	80	97.5	+	0	0	I13MI24MD42M	s0;size=5	s1
H	9	80	93.8	+	0	0	5D75M12I	s0;size=5	s9
H	7	81	97.5	+	0	0	D13MD24MI42M	s1	s7
H	6	81	96.2	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	92.4	+	0	0	D13MD66M	s1	s8
H	9	81	91.4	+	0	0	6D8MD24MI42M12I	s1	s9
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	97.5	+	0	0	40MD39M	s6	s8
H	9	80	92.5	+	0	0	5D75M12I	s6	s9
H	8	80	98.7	+	0	0	41MD38M	s7	s8
H	9	80	93.8	+	0	0	5D75M12I	s7	s9
H	9	79	92.4	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn takes_identity_definition_1() {
    let args = ["--id", "0.9", "--iddef", "1"];
    assert_eq!(
        allpr(&args, &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	95.1	+	0	0	I13MI24MD42M	s0;size=5	s1
H	7	81	95.1	+	0	0	D13MD24MI42M	s1	s7
H	6	81	93.9	+	0	0	D13MD24MI42M	s1	s6
H	8	81	90.1	+	0	0	D13MD66M	s1	s8
N	*	*	*	.	*	*	*	s2	*
N	*	*	*	.	*	*	*	s3;size=2	*
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	8	80	97.5	+	0	0	41MD38M	s7	s8
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn takes_identity_definition_3() {
    let args = ["--id", "0.9", "--iddef", "3"];
    assert_eq!(
        allpr(&args, &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	9	80	97.7	+	0	0	5D75M12I	s0;size=5	s9
H	3	80	97.5	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	95.1	+	0	0	I13MI24MD42M	s0;size=5	s1
H	9	81	95.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	95.1	+	0	0	D13MD24MI42M	s1	s7
H	6	81	93.8	+	0	0	D13MD24MI42M	s1	s6
H	3	81	92.6	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	90.1	+	0	0	D13MD66M	s1	s8
N	*	*	*	.	*	*	*	s2	*
H	9	60	97.7	+	0	0	5I60M22I	s3;size=2	s9
H	7	60	97.5	+	0	0	10I60M10I	s3;size=2	s7
H	6	60	96.2	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	94.9	+	0	0	10I31MD28M10I	s3;size=2	s8
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	96.6	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	97.7	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	95.4	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn takes_identity_definition_4() {
    let args = ["--id", "0.9", "--iddef", "4"];
    assert_eq!(
        allpr(&args, &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	95.1	+	0	0	I13MI24MD42M	s0;size=5	s1
H	7	81	95.1	+	0	0	D13MD24MI42M	s1	s7
H	6	81	93.9	+	0	0	D13MD24MI42M	s1	s6
H	8	81	90.1	+	0	0	D13MD66M	s1	s8
N	*	*	*	.	*	*	*	s2	*
N	*	*	*	.	*	*	*	s3;size=2	*
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	8	80	97.5	+	0	0	41MD38M	s7	s8
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn keeps_the_first_hits() {
    let main = inputs().main;
    let args = ["--accept-all", "--maxhits", "2"];
    assert_eq!(
        allpr(&args, &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
    let args = ["--id", "0.5", "--top-hits-only"];
    assert_eq!(
        allpr(&args, &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn writes_no_weak_hits() {
    let args = ["--id", "0.95", "--weak-id", "0.8"];
    assert_eq!(
        allpr(&args, &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn counts_n_as_a_mismatch() {
    let args = ["--id", "0.9", "--n-mismatch"];
    assert_eq!(
        allpr(&args, &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	6	80	96.2	+	0	0	80M	s0;size=5	s6
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	6	81	92.6	+	0	0	D13MD24MI42M	s1	s6
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	6	60	95.0	+	0	0	10I60M10I	s3;size=2	s6
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	96.2	+	0	0	80M	s6	s7
H	9	80	96.0	+	0	0	5D75M12I	s6	s9
H	8	80	95.0	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxdiffs() {
    assert_eq!(
        allpr(&["--id", "0", "--max-diffs", "2"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxsubs() {
    assert_eq!(
        allpr(&["--id", "0", "--max-subs", "0"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
N	*	*	*	.	*	*	*	s6	*
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxgaps() {
    assert_eq!(
        allpr(&["--id", "0", "--max-gaps", "1"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_mincols() {
    assert_eq!(
        allpr(&["--id", "0", "--min-cols", "70"], &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
N	*	*	*	.	*	*	*	s3;size=2	*
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_leftjust() {
    assert_eq!(
        allpr(&["--id", "0", "--left-just"], &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
N	*	*	*	.	*	*	*	s3;size=2	*
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	8	80	97.5	+	0	0	41MD38M	s7	s8
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_rightjust() {
    assert_eq!(
        allpr(&["--id", "0", "--right-just"], &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
N	*	*	*	.	*	*	*	s3;size=2	*
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	8	80	97.5	+	0	0	41MD38M	s7	s8
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_qcov() {
    assert_eq!(
        allpr(&["--id", "0", "--query-cov", "0.9"], &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_tcov() {
    assert_eq!(
        allpr(&["--id", "0", "--target-cov", "0.9"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
N	*	*	*	.	*	*	*	s3;size=2	*
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	8	80	97.5	+	0	0	41MD38M	s7	s8
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxid() {
    assert_eq!(
        allpr(&["--id", "0", "--max-id", "0.97"], &inputs().main),
        r##"H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
N	*	*	*	.	*	*	*	s7	*
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_mid() {
    assert_eq!(
        allpr(&["--id", "0", "--mid", "99"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
N	*	*	*	.	*	*	*	s6	*
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_idprefix() {
    assert_eq!(
        allpr(&["--id", "0", "--id-prefix", "5"], &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
N	*	*	*	.	*	*	*	s1	*
N	*	*	*	.	*	*	*	s2	*
N	*	*	*	.	*	*	*	s3;size=2	*
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	8	80	97.5	+	0	0	41MD38M	s7	s8
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_idsuffix() {
    assert_eq!(
        allpr(&["--id", "0", "--id-suffix", "3"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
H	7	80	98.8	+	0	0	80M	s6	s7
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	8	80	97.5	+	0	0	41MD38M	s7	s8
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxqt() {
    assert_eq!(
        allpr(&["--id", "0", "--max-qt", "0.95"], &inputs().main),
        r##"H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_minqt() {
    assert_eq!(
        allpr(&["--id", "0", "--min-qt", "1.05"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
N	*	*	*	.	*	*	*	s3;size=2	*
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
N	*	*	*	.	*	*	*	s6	*
N	*	*	*	.	*	*	*	s7	*
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_minsl() {
    assert_eq!(
        allpr(&["--id", "0", "--min-sl", "0.9"], &inputs().main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
N	*	*	*	.	*	*	*	s3;size=2	*
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxsl() {
    assert_eq!(
        allpr(&["--id", "0", "--max-sl", "0.95"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxqsize() {
    assert_eq!(
        allpr(&["--id", "0", "--max-qsize", "2"], &inputs().main),
        r##"N	*	*	*	.	*	*	*	s0;size=5	*
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_mintsize() {
    assert_eq!(
        allpr(&["--id", "0", "--min-tsize", "2"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
N	*	*	*	.	*	*	*	s3;size=2	*
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
N	*	*	*	.	*	*	*	s6	*
N	*	*	*	.	*	*	*	s7	*
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_minsizer() {
    assert_eq!(
        allpr(&["--id", "0", "--min-size-ratio", "2"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
N	*	*	*	.	*	*	*	s1	*
N	*	*	*	.	*	*	*	s2	*
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
N	*	*	*	.	*	*	*	s4	*
N	*	*	*	.	*	*	*	s5	*
N	*	*	*	.	*	*	*	s6	*
N	*	*	*	.	*	*	*	s7	*
N	*	*	*	.	*	*	*	s8	*
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_maxsizer() {
    assert_eq!(
        allpr(&["--id", "0", "--max-size-ratio", "1"], &inputs().main),
        r##"N	*	*	*	.	*	*	*	s0;size=5	*
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
N	*	*	*	.	*	*	*	s3;size=2	*
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn applies_selfid() {
    assert_eq!(
        allpr(&["--id", "0", "--self-id"], &inputs().main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn skips_targets_of_the_same_label() {
    let selves = inputs().selves;
    assert_eq!(
        allpr(&["--id", "0", "--self"], &selves),
        r##"H	2	80	98.8	+	0	0	50MD29M	a	b
H	2	82	96.3	+	0	0	4M2D46MD29M	a	b
N	*	*	*	.	*	*	*	b	*
"##
    );
    assert_eq!(
        allpr(&["--id", "0"], &selves),
        r##"H	2	80	98.8	+	0	0	50MD29M	a	b
H	1	80	97.6	+	0	0	4M2I76M	a	a
H	2	82	96.3	+	0	0	4M2D46MD29M	a	b
N	*	*	*	.	*	*	*	b	*
"##
    );
}

#[test]
fn takes_the_scores() {
    let main = inputs().main;
    let args = ["--accept-all", "--match", "1", "--mismatch", "-2"];
    assert_eq!(
        allpr(&args, &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	37.7	+	0	0	3D77M3I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	81.2	+	0	0	2D37MI26MD15M	s1	s2
H	4	81	54.4	+	0	0	2D12MD66M2I	s1	s4
H	5	81	36.1	+	0	0	9D72M8I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	52.6	+	0	0	D78M2I	s2	s4
H	5	79	38.4	+	0	0	6D73M7I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	38.3	+	0	0	7I60M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	55.1	+	0	0	I78M2D	s4	s8
H	5	80	33.8	+	0	0	80M	s4	s5
H	9	80	38.5	+	0	0	2D78M9I	s5	s9
H	6	80	37.7	+	0	0	3I77M3D	s5	s6
H	7	80	37.7	+	0	0	3I77M3D	s5	s7
H	8	80	32.9	+	0	0	3I76M4D	s5	s8
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
    let args = ["--accept-all", "--gapopen", "10I/1E", "--gapext", "1I/1E"];
    assert_eq!(
        allpr(&args, &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	90.1	+	0	0	2MI6MD55MD15M	s0;size=5	s2
H	4	80	79.3	+	0	0	4MD25MI20MI11M3D16M2I	s0;size=5	s4
H	5	80	36.6	+	0	0	4D2M10I5M3D4M27I11MD6M7D5M2D10M20D	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	3	81	96.7	+	0	0	11D3MD24MI32M10D	s1	s3;size=2
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	8	81	93.8	+	0	0	D13MD24MI3MD38M	s1	s8
H	2	81	87.8	+	0	0	9MD2MD26MI26MD15M	s1	s2
H	4	81	75.9	+	0	0	D4MD8MD17MI7MI12MI11M3D16M2I	s1	s4
H	5	81	36.6	+	0	0	5D2M10I3M4D6M27I11MD6M6D5M2D10M20D	s1	s5
H	6	79	91.4	+	0	0	2MD6MI55MI15M	s2	s6
H	7	79	90.1	+	0	0	2MD6MI55MI15M	s2	s7
H	8	79	90.1	+	0	0	2MD6MI31MD23MI15M	s2	s8
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	3	79	79.7	+	0	0	10D56M9D4M	s2	s3;size=2
H	4	79	71.6	+	0	0	2MD27MI20MI11M2D16M2I	s2	s4
H	5	79	41.1	+	0	0	6D2M9I10M9I2M15D17M3D4M4I20M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	76.1	+	0	0	9I20MI20MI11MI3M10I4M2D	s3;size=2	s4
H	5	60	37.7	+	0	0	13I8M27I11MD6M7D5M2D10M10D	s3;size=2	s5
H	7	80	79.3	+	0	0	4MI25MD20MD11M3I16M2D	s4	s7
H	6	80	78.0	+	0	0	4MI25MD20MD11M3I16M2D	s4	s6
H	8	80	76.8	+	0	0	4MI25MD11MD7MD12M3I16M2D	s4	s8
H	9	80	70.0	+	0	0	4D25MD20MD11MD3M16I12M2D	s4	s9
H	5	80	33.0	+	0	0	7I10M3D4M25I19M17D9M10D6M2D	s4	s5
H	8	80	37.0	+	0	0	4I2M10D5M3I4M27D11MI6M6I4M2I11M20I	s5	s8
H	7	80	36.6	+	0	0	4I2M10D5M3I4M27D11MI6M7I5M2I10M20I	s5	s7
H	6	80	35.9	+	0	0	21D7M18D15M17I9M2I10M20I	s5	s6
H	9	80	30.3	+	0	0	4I7M3D3M10I6M13D7M35I15M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
    let args = ["--accept-all", "--gapopen", "30Q/5T", "--gapext", "3QI/1TI"];
    assert_eq!(
        allpr(&args, &main),
        r##"H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	MI12MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	3	80	85.7	+	0	0	10D56M10D4M	s0;size=5	s3;size=2
H	4	80	76.5	+	0	0	4MD25MI19M3D4M4DMD2M4D2MD6M2DM15I	s0;size=5	s4
H	5	80	44.7	+	0	0	4D5M3D15M2DMD6MDM8D6M4D4M6D13M29I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	MD12MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	MD12MD24MI42M	s1	s6
H	8	81	91.4	+	0	0	MD12MD66M	s1	s8
H	2	81	87.8	+	0	0	9MD2MD26MI26MD15M	s1	s2
H	3	81	83.1	+	0	0	11D3MD24MI28M10D4M	s1	s3;size=2
H	4	81	60.5	+	0	0	MD3MD8MD17M10D4M2D2M2D2M4D2M4DM5D6M2D3M31I	s1	s4
H	5	81	42.0	+	0	0	M9D19M2DMD7M8D6M4D4M6D13M29I	s1	s5
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	9	79	86.2	+	0	0	MD2M2D5M2D51MI15M12I	s2	s9
H	3	79	75.0	+	0	0	3D2M5DM2D53M9D4M	s2	s3;size=2
H	4	79	68.8	+	0	0	2MD27MI19M3D4M4DMD2M3D2MD6M2DM15I	s2	s4
H	5	79	44.3	+	0	0	MD3M6D16M2DMD6MDM8D6M4D4M5D13M29I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	78.7	+	0	0	9I20MI19M3D4M4DMD2M4D2M22I	s3;size=2	s4
H	5	60	50.0	+	0	0	47I4DM3D11MD6M7D5M2D10M10D	s3;size=2	s5
H	7	80	77.8	+	0	0	4MI25MD20MD11MD3M10D4M12I	s4	s7
H	9	80	77.6	+	0	0	4D25MD20MD11MD3M10D4M24I	s4	s9
H	6	80	76.5	+	0	0	4MI25MD20MD11MD3M10D4M12I	s4	s6
H	8	80	75.3	+	0	0	4MI25MD11MD7MD12MD3M10D4M12I	s4	s8
H	5	80	46.2	+	0	0	21I6M2DM4D7M4D2MD15M7D2M2D9MD17M	s4	s5
H	7	80	43.8	+	0	0	M8D4MD15MD22M6D5M2D2M2D11M20I	s5	s7
H	8	80	43.8	+	0	0	M8D4MD15MD19M7D7M2D2M2D11M20I	s5	s8
H	9	80	43.8	+	0	0	2M8D11M6D2MD22M6D5M2D2M2D11M32I	s5	s9
H	6	80	40.0	+	0	0	35I2M12D9M4D7MD2M6D4MD8M6DM2D6M3D6M	s5	s6
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn forbids_infinite_gaps() {
    // the 16-bit aligner hands every pair to the linear-memory one; with
    // --id, hits that open an infinite gap are rejected
    let main = inputs().main;
    let args = ["--accept-all", "--gapopen", "*I/2E"];
    assert_eq!(
        allpr(&args, &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	2	80	77.2	+	0	0	79MD	s0;size=5	s2
H	1	80	76.2	+	0	0	I80M	s0;size=5	s1
H	8	80	69.6	+	0	0	79MD	s0;size=5	s8
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	4	80	50.6	+	0	0	D79MI	s0;size=5	s4
H	8	81	81.0	+	0	0	2D79M	s1	s8
H	7	81	76.2	+	0	0	D80M	s1	s7
H	9	81	76.0	+	0	0	6D75M12I	s1	s9
H	6	81	75.0	+	0	0	D80M	s1	s6
H	3	81	70.0	+	0	0	11D60M10D	s1	s3;size=2
H	2	81	60.8	+	0	0	2D79M	s1	s2
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	4	81	44.9	+	0	0	3D78M2I	s1	s4
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	9	79	79.7	+	0	0	5D74M13I	s2	s9
H	6	79	78.5	+	0	0	79MI	s2	s6
H	7	79	77.2	+	0	0	79MI	s2	s7
H	8	79	67.1	+	0	0	79M	s2	s8
H	5	79	58.8	+	0	0	62D17M63I	s2	s5
H	4	79	52.6	+	0	0	D78M2I	s2	s4
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	68.3	+	0	0	10I60M9I	s3;size=2	s8
H	4	60	51.7	+	0	0	9I60M11I	s3;size=2	s4
H	5	60	38.3	+	0	0	7I60M13I	s3;size=2	s5
H	8	80	55.1	+	0	0	I78M2D	s4	s8
H	9	80	52.6	+	0	0	4D76M11I	s4	s9
H	6	80	51.9	+	0	0	I79MD	s4	s6
H	7	80	50.6	+	0	0	I79MD	s4	s7
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	8	80	40.7	+	0	0	21D59M20I	s5	s8
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	70.9	+	0	0	79MD	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	69.6	+	0	0	79MD	s7	s8
H	9	79	67.6	+	0	0	5D74M13I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
    let args = ["--id", "0", "--gapopen", "*I/2E"];
    assert_eq!(
        allpr(&args, &main),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	2	80	77.2	+	0	0	79MD	s0;size=5	s2
H	1	80	76.2	+	0	0	I80M	s0;size=5	s1
H	8	80	69.6	+	0	0	79MD	s0;size=5	s8
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	4	80	50.6	+	0	0	D79MI	s0;size=5	s4
H	8	81	81.0	+	0	0	2D79M	s1	s8
H	7	81	76.2	+	0	0	D80M	s1	s7
H	9	81	76.0	+	0	0	6D75M12I	s1	s9
H	6	81	75.0	+	0	0	D80M	s1	s6
H	3	81	70.0	+	0	0	11D60M10D	s1	s3;size=2
H	2	81	60.8	+	0	0	2D79M	s1	s2
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	4	81	44.9	+	0	0	3D78M2I	s1	s4
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	9	79	79.7	+	0	0	5D74M13I	s2	s9
H	6	79	78.5	+	0	0	79MI	s2	s6
H	7	79	77.2	+	0	0	79MI	s2	s7
H	8	79	67.1	+	0	0	79M	s2	s8
H	5	79	58.8	+	0	0	62D17M63I	s2	s5
H	4	79	52.6	+	0	0	D78M2I	s2	s4
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	68.3	+	0	0	10I60M9I	s3;size=2	s8
H	4	60	51.7	+	0	0	9I60M11I	s3;size=2	s4
H	5	60	38.3	+	0	0	7I60M13I	s3;size=2	s5
H	8	80	55.1	+	0	0	I78M2D	s4	s8
H	9	80	52.6	+	0	0	4D76M11I	s4	s9
H	6	80	51.9	+	0	0	I79MD	s4	s6
H	7	80	50.6	+	0	0	I79MD	s4	s7
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	8	80	40.7	+	0	0	21D59M20I	s5	s8
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	70.9	+	0	0	79MD	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	69.6	+	0	0	79MD	s7	s8
H	9	79	67.6	+	0	0	5D74M13I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn aligns_long_sequences() {
    let long = inputs().long;
    assert_eq!(
        allpr(&["--accept-all"], &long),
        r##"H	1	4000	99.3	+	0	0	173MD68MI115MI172MI349MD169MI231MI253MI137MD251MI162MD490MD24MI106MI206MD697MD390M	l0	l1
H	2	4000	40.6	+	0	0	132D22M3D90M8D39M6I66MD111M3D65M3D42M9I35M10I35M3D6M3D27M5D28MD63M4D27M3I12M2D68M11I89M35I36M10D10M8D30M9I16M6I53MD36MD76M16D38M5I122M5I13M3D47M13D26M39D47MD37M4I50M3I30MI26M15D58M2I27M8I139MD59M4D39M9I36M2D25M8D31M7D9M2I91M3I15M9I19M4D31MD50M8I131M14D17M19D35M18D24M6D17MD19M10D14MD36M2D17M2I13M7D10MD28M3D10M2I31M7D8M2D43M4I100M3D26M3I8M4I10MI8M4I36M9I130M19I36M17I83MD65M6I26M13D42M7I23MD33M7I59MD96M6D43M6D42M12I31M22D57M191I	l0	l2
H	2	4002	41.3	+	0	0	44D13MD16M3I30M10I18MI25M4I22M7D35M12D71M9D41M8D23M6D74M24I79MD24M3D27M5D44M5D16M3I23M3I15MD27MD91M13I50M5D9M4D44M4D73M5I14MI26M2I24M10I77M4I29M9I58MD22MI26M28D67M17D6M3D43M4D15M6D24M3D35M4D13M12D24M2D5M5D42M3D121MD11M2I19M4D38M3D58MI57M8I122M7I13M4I27M8D24MI18M4I52M17I18MD30M2D42M11I35MD28M3I81M15D19MD17M5I41M3D68MI42M3I21M2I41M9I19MI36MD28M7D17M12D53M5D18M7I59MD116M22I19M10I32M7I37M5I28M10I22MI46M14I55M9D21M5D26M5D60M6I15M4I42M2I32M11D68M2I54M2D74M30D37M7D36M11I40MD17M9I35MI38M44I	l1	l2
N	*	*	*	.	*	*	*	l2	*
"##
    );
    // scores beyond 16 bits, high or low, fall back to the linear-memory
    // aligner
    let args = ["--accept-all", "--match", "10"];
    assert_eq!(
        allpr(&args, &long),
        r##"H	1	4000	99.4	+	0	0	173MD68MI116MI171MI349MD169MI231MI253MI137MD252MI161MD490MD24MI106MI206MD662MI8MD26MD390M	l0	l1
H	2	4000	48.1	+	0	0	11D5MI6M3I5MI14M10D6MD6MI6MD7M6D9M7I11M8I7M8I5M3I5MD9MD15M3D15M2D7M7I16M2D2M2D9M6I19M3I6MD10MD4MD14M10D14M3I26MD10M2D21MD11MD6M4D9M3D3M5D25MI17M5D4MD7M3I17MI10MI21M6D5M17D7MI12M3D12M3D13M6D5M3D6M3I2MI8M2D6M2I12M2D10M3D11M2D2MI4M4I23M3I15MD27MD66M3D13MD5M2I10MI10M2D6M3D15M2D6M9D12M4D17M3D10M9D4M2D13MD14M5D8M2I14M2D36M2D6MI7M3D13M3I10M7D9MD16MD5MD21M2I13M3D22M2I6M3I10MI5M2I3MI6M3I6M5D50M2D3M2D6M7D7M3D3M2D7MD17M4I4M2I32M18I12M2I5M8I12MI6MI8M6I4M4I13MI24M2I12M2D8MD3MD2M3D12M8I6M8I9M6I14MI3M4I4M3I10M3I7M6I9M5D11MI2M2D14M8I17M5D9M6I8MI12MI4M5I4MD20MI6M5D18M3I12M2I13M8I4M3I16M14I7MD6M3I5MI32M4I15M2D10MI3M3D12M8I7M3D12MI7M3D27MD5M5D5MD6MI16M2D15M6I4M2D7M2I20M2D6M2I3MD14M4I12MD11M2D4M16D6M3D9M6D9M12D8M15D10M4D15M3D12M4D12M6D6M5D12MD9MD3M3D7M6D27M4D9M2I29M9I7M7I8MD3M3D9M3I11MD11M14I17M4D3MD10M7D15MD4M13D22MD7MI12M4D8MD17M5I41MD4M5D8M3I5MD5MI2M4I10MD10M7I9M3I11M2I3M6I7M3I9M2I22M7D14M5I4MD7M4D8M8I4M5I15M8D10M3I10MD8M2D14MD8M7D12M7D17M7D6MI4M6D42M8I10MD8M4D21M5I5M6I4MI11M3I15M10I2MD12M4I7M3D8M2I6M2I24M3D7MI9M16I3M4I3M2I5MD3M2D13M4I15MI4M2I10M4I10MD13MI9M3I14M2D27M3D11M3D2MD22MI31MD10M2I25MI2MI9M2D7M2D5MD6MI3M2I3M3I3MI21M5D26M5I3M4I5MI8MI4M4I4MD11M4I5M3I3MD8M2D6M7D5MD7M4D12MD24M2I32M2D19M9D23MD26M2I23M4D7M3I11MI9M2D26M6D4M4D3M2I12M7I11M7D20M3I10M4I13M4D6MI15MI9M2I22MI8M2I8M3I4MD5MD11MI19M8I15M2D17M2D15MI8M5I27M13I	l0	l2
H	2	4002	48.3	+	0	0	11D5MI6M3I5MI14M10D6MD6MI6MD7M6D9M7I11M8I7M8I5M3I5MD9MD15M4I18MD23M2D2M2D9M6I19M5I6M5D12MD14M10D14M3I26MD10M2D21MD11M2D6M4D9M3D3M5D25MI17M5D4MD7M3I17MI10MI21M6D5M17D20M3D12M3D13M6D5M3D6M3I2MI8M2D6M2I12M2D10M3D11M2D2MI4M4I23M3I15MD27MD66M3D13MD5M2I10MI10M11I2M3I31M5D9M4D37MD9M3D71M4D2M5I5M3I4M2I25M2I16M7I14M2D11M4I7MI51M4I7M2D7M2I8M2I2M4I6M3I7M3I33M3D12MD10M5D14M11D7MI9M3D9M9D29M2I7M2D30MI3M3I37MD10M2D5M2D17MI8M2I3M3I4M10I9M5I15M12I22M5I7M3D20M7I6M7I9MI4M3I13MD6M4I6MI20M5D7MI9M3D14M9D24M3D7M2I18M4D21M14I7MD6M3I5MI32M4I15M2D10MI3M3D12M8I7M3D12MI7M3D27MD5M5D5MD6MI5M2D11MD15M6I4M2D7M2I20M2D6M2I3MD14M4I12MD11M2D4M6I6M2D8MI5MI5M2D8M6I3M11I19M6I21M5I10MD13M5I2M2I11MD16M2I3M4D14M3D9M5D7M3D10MD16MI12M6I21M2I2MI14MD3MD43MI3MI11M2I13MI4MI6M6I18M2I10MD3MI12M3I4M2D18M6D9M2D9MD11M5I15M6I14M4I13M10I3MI6M3D13M2D18M3I6M2D7M2D13M2D22M4D12M4D3MI13M2I14M5I19M2D11M6D7M12D6M3D4M5D7M2I36M2D16M7D18M3D15MD10MD4MI2M7D14M10D4M2I4M2I2M2I4MD7M4D14M9D12M5D26MD25M2D6M6D7M4D23MD9M9D14MI4M2I10M4I10MD13MI9M3I14M2D27M3D11M3D2MD22MI31MD10M2I25MI2MI9M2D7M2D5MD6MI3M2I3M3I3MI21M5D26M5I3M4I5MI8MI4M4I4MD11M4I5M3I3MD8M2D6M7D5MD7M4D12MD24M2I32M2D19M4D12M5D37M2I23M4D7M3I11MI9M2D26M6D4M4D3M2I12M7I11M7D20M3I10M4I13M4D6MI15MI9M2I22MI8M2I8M3I4MD5MD11MI19M8I15M2D17M2D15MI8M5I27M13I	l1	l2
N	*	*	*	.	*	*	*	l2	*
"##
    );
    let args = ["--accept-all", "--mismatch", "-1000"];
    assert_eq!(
        allpr(&args, &long),
        r##"H	1	4000	99.2	+	0	0	173MD68MI115MI172MI299MID49MD63MID72MI2MD30MI231MI134MID18MID99MI137MD243MID7MI162MD490MD24MI106MI76MID129MD662MI8MD26MD253MID136M	l0	l1
H	2	4000	0.0	+	0	0	3990I4000D	l0	l2
H	2	4002	0.0	+	0	0	3990I4002D	l1	l2
N	*	*	*	.	*	*	*	l2	*
"##
    );
    // so do sequences too long for it
    assert_eq!(
        allpr(&["--accept-all"], &inputs().big),
        r##"H	1	5100	99.6	+	0	0	182MI358MI9MI91MD231MI502MI159MI319MI383MI932MI441MI51MI183MI238MD136MD527MD52MD301M	g0	g1
N	*	*	*	.	*	*	*	g1	*
"##
    );
}

#[test]
fn aligns_empty_sequences() {
    // vsearch with --minseqlength 0
    let input = ">e0\n\n>e1\nACGTACGT\n>e2\n\n>e3\nACGAACGT\n";
    assert_eq!(
        allpr(&["--accept-all"], input),
        r##"H	1	0	0.0	+	0	0	8I	e0	e1
H	2	0	0.0	+	0	0	=	e0	e2
H	3	0	0.0	+	0	0	8I	e0	e3
H	3	8	87.5	+	0	0	8M	e1	e3
H	2	8	0.0	+	0	0	8D	e1	e2
H	3	0	0.0	+	0	0	8I	e2	e3
N	*	*	*	.	*	*	*	e3	*
"##
    );
}

#[test]
fn reads_fastq() {
    let fastq: String = inputs()
        .main
        .split('>')
        .skip(1)
        .map(|record| {
            let (header, seq) = record.split_once('\n').unwrap();
            let seq = seq.trim_end();
            format!("@{header}\n{seq}\n+\n{}\n", "I".repeat(seq.len()))
        })
        .collect();
    assert_eq!(
        allpr(&["--accept-all"], &fastq),
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn checks_abundances() {
    // every sequence is read before any is aligned
    let output = run(&["--accept-all"], ">a\nACGT\n>b;size=0\nACGT\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("abundance"), "{stderr}");
}

#[test]
fn rejects_invalid_arguments() {
    let main = inputs().main;
    for args in [
        &[][..],
        &["--id", "1.5"],
        &["--id", "-0.1"],
        &["--id", "0.9", "--accept-all"],
        &["--id", "0.9", "--iddef", "5"],
        &["--id", "0.9", "--maxhits", "0"],
        &["--id", "0.9", "--gapopen", "2EL"],
    ] {
        assert_eq!(run(args, &main).status.code(), Some(2), "{args:?}");
    }
}

#[test]
fn reads_and_writes_files() {
    let dir = std::env::temp_dir();
    let input = dir.join(format!("rsearch-allpr-{}-in.fa", std::process::id()));
    let out = dir.join(format!("rsearch-allpr-{}-out.uc", std::process::id()));
    std::fs::write(&input, inputs().main).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_allpr"))
        .args([
            "--accept-all",
            "-i",
            input.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&out);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(
        written,
        r##"H	3	80	100.0	+	0	0	10D60M10D	s0;size=5	s3;size=2
H	7	80	100.0	+	0	0	=	s0;size=5	s7
H	9	80	100.0	+	0	0	5D75M12I	s0;size=5	s9
H	6	80	98.8	+	0	0	80M	s0;size=5	s6
H	8	80	97.5	+	0	0	41MD38M	s0;size=5	s8
H	1	80	96.3	+	0	0	I13MI24MD42M	s0;size=5	s1
H	2	80	86.2	+	0	0	64MD15M	s0;size=5	s2
H	4	80	75.0	+	0	0	D29MI20MI29MD	s0;size=5	s4
H	5	80	58.8	+	0	0	63D17M63I	s0;size=5	s5
H	9	81	97.4	+	0	0	6D8MD24MI42M12I	s1	s9
H	7	81	96.3	+	0	0	D13MD24MI42M	s1	s7
H	6	81	95.1	+	0	0	D13MD24MI42M	s1	s6
H	3	81	93.3	+	0	0	12D27MI32M10D	s1	s3;size=2
H	8	81	91.2	+	0	0	D13MD66M	s1	s8
H	2	81	85.4	+	0	0	9M2D28MI26MD15M	s1	s2
H	4	81	72.8	+	0	0	2D12MD17MI7MI12MI29MD	s1	s4
H	5	81	58.8	+	0	0	64D17M63I	s1	s5
H	9	79	89.3	+	0	0	5D59MI15M12I	s2	s9
H	6	79	87.5	+	0	0	64MI15M	s2	s6
H	3	79	86.7	+	0	0	10D60M9D	s2	s3;size=2
H	7	79	86.2	+	0	0	64MI15M	s2	s7
H	8	79	86.2	+	0	0	40MD23MI15M	s2	s8
H	4	79	70.0	+	0	0	D29MI20MI11M2D16M2I	s2	s4
H	5	79	48.1	+	0	0	7D19M4I28MI25M3I	s2	s5
H	7	60	100.0	+	0	0	10I60M10I	s3;size=2	s7
H	9	60	100.0	+	0	0	5I60M22I	s3;size=2	s9
H	6	60	98.3	+	0	0	10I60M10I	s3;size=2	s6
H	8	60	96.7	+	0	0	10I31MD28M10I	s3;size=2	s8
H	4	60	82.3	+	0	0	9I20MI20MI20M9I	s3;size=2	s4
H	5	60	42.2	+	0	0	3I16M4I44M13I	s3;size=2	s5
H	9	80	77.6	+	0	0	4D25MD20MD29M13I	s4	s9
H	7	80	75.0	+	0	0	I29MD20MD29MI	s4	s7
H	6	80	73.8	+	0	0	I29MD20MD29MI	s4	s6
H	8	80	68.8	+	0	0	I29MD13M2D35MI	s4	s8
H	5	80	45.0	+	0	0	40I40M40D	s4	s5
H	6	80	58.8	+	0	0	63I17M63D	s5	s6
H	7	80	58.8	+	0	0	63I17M63D	s5	s7
H	8	80	58.8	+	0	0	62I17M63D	s5	s8
H	9	80	42.6	+	0	0	33I54M26D	s5	s9
H	7	80	98.8	+	0	0	80M	s6	s7
H	9	80	98.7	+	0	0	5D75M12I	s6	s9
H	8	80	96.2	+	0	0	40MD39M	s6	s8
H	9	80	100.0	+	0	0	5D75M12I	s7	s9
H	8	80	97.5	+	0	0	41MD38M	s7	s8
H	9	79	97.3	+	0	0	5D36MI38M12I	s8	s9
N	*	*	*	.	*	*	*	s9	*
"##
    );
}

#[test]
fn writes_the_same_hits_on_any_thread_count() {
    let input = many_sequences(400);
    for args in [
        &["--id", "0.8"][..],
        &["--accept-all", "--maxhits", "3"],
        &["--id", "0.7", "--top-hits-only", "--iddef", "1"],
        &["--id", "0.9", "--max-diffs", "4", "--query-cov", "0.9"],
    ] {
        let one = same_on_threads(args, &input, &["2", "3", "5", "16", "0"]);
        assert!(one.status.success(), "{args:?}");
        assert!(!one.stdout.is_empty(), "{args:?}");
    }
}

#[test]
fn stops_at_an_invalid_sequence_as_on_one_thread() {
    let input = many_sequences(300);
    for bad in [">bad;size=0\nACGT\n", ">bad\nAC\u{1}GT\n"] {
        for whole in [format!("{input}{bad}{input}"), format!("{bad}{input}")] {
            let one = same_on_threads(&["--id", "0.8"], &whole, &["2", "3", "8"]);
            assert_eq!(one.status.code(), Some(1));
        }
    }
}

#[test]
fn takes_the_thread_count() {
    let main = inputs().main;
    let expected = allpr(&["--id", "0.8"], &main);
    for t in [&["-t", "0"][..], &["--threads", "5"], &["-t", "1025"]] {
        let args: Vec<&str> = ["--id", "0.8"]
            .iter()
            .copied()
            .chain(t.iter().copied())
            .collect();
        let output = run(&args, &main);
        assert!(output.status.success(), "{t:?}");
        assert_eq!(String::from_utf8(output.stdout).unwrap(), expected, "{t:?}");
    }
    for t in ["-1", "x", "1.5", ""] {
        assert_eq!(
            run(&["--id", "0.8", "-t", t], &main).status.code(),
            Some(2),
            "{t:?}"
        );
    }
    assert_eq!(
        same_on_threads(&["--id", "0.8"], "", &["2", "8"]).stdout,
        b""
    );
}

/// Runs allpr with `args` on `input`, fed on stdin.
fn run(args: &[&str], input: &str) -> Output {
    run_at(None, args, input)
}

/// The SIMD levels of this processor, scalar with SSE2 and AVX2 or with
/// NEON. Each aligns eight or sixteen targets at once and must give the
/// output of the others.
fn levels() -> Vec<rsearch::simd::Level> {
    let levels = rsearch::simd::available();
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    assert!(levels.len() >= 2, "{levels:?}");
    levels
}

#[test]
fn writes_the_same_hits_at_every_simd_level() {
    let levels = levels();
    let input = many_sequences(400);
    for args in [
        &["--id", "0.8"][..],
        &["--accept-all", "--maxhits", "3"],
        &["--id", "0.7", "--top-hits-only", "--iddef", "1"],
    ] {
        let runs: Vec<Output> = levels
            .iter()
            .map(|&level| run_at(Some(level), args, &input))
            .collect();
        assert!(runs[0].status.success(), "{args:?}");
        assert!(!runs[0].stdout.is_empty(), "{args:?}");
        for (level, output) in levels.iter().zip(&runs).skip(1) {
            assert!(output.status.success(), "{level:?} {args:?}");
            assert!(
                output.stdout == runs[0].stdout,
                "{level:?} writes other hits, {args:?}"
            );
        }
    }
}
