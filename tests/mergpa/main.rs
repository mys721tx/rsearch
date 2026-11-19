//! End-to-end tests of the mergpa binary.
//!
//! Expected outputs are those of `vsearch --fastq_mergepairs` with
//! `--reverse` on the same inputs. Merged pairs go to `-o` (`--fastqout`).
//! The forward and reverse reads of pairs not merged go to `--not-merged`
//! and `--rev-not-merged` (`--fastqout_notmerged_fwd`,
//! `--fastqout_notmerged_rev`). vsearch's merging report goes to `--report`
//! (what vsearch prints on stderr), and the errors of each merged pair to
//! `--ee-table` (`--eetabbedout`). `--min-overlap`, `--max-diffs`,
//! `--max-diff-pct` and `--stagger` are `--fastq_minovlen`,
//! `--fastq_maxdiffs`, `--fastq_maxdiffpct` and
//! `--fastq_allowmergestagger`. `--ascii`, `--qmin`, `--qmax`, `--qmin-out`
//! and `--qmax-out` are `--fastq_ascii`, `--fastq_qmin`, `--fastq_qmax`,
//! `--fastq_qminout` and `--fastq_qmaxout`. Truncation and the read filters
//! (`--fastq_truncqual`, `--fastq_minlen`, `--fastq_maxlen`,
//! `--fastq_maxns`) are `sqtrim` and `sqfilt` before mergpa, and the merged
//! filters (`--fastq_maxee`, `--fastq_minmergelen`, `--fastq_maxmergelen`)
//! `sqfilt` after it.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-mergpa-{}-{name}", std::process::id()))
}

/// Builds FASTQ text from (label, sequence, quality) triples.
fn fastq(records: &[(&str, &str, &str)]) -> String {
    records
        .iter()
        .map(|(id, seq, qual)| format!("@{id}\n{seq}\n+\n{qual}\n"))
        .collect()
}

/// The reverse complement of an upper-case sequence.
fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

/// A fragment of 80 bases: pairs of reads of 50 from each end overlap by 20.
const FRAGMENT: &str =
    "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGT";

/// A fragment of 104 bases: pairs of reads of 70 overlap by 36.
fn long_fragment() -> String {
    format!("{FRAGMENT}CCATGACTTGGAACCTGTAGACTC")
}

/// `n` copies of a quality symbol.
fn q(symbol: char, n: usize) -> String {
    symbol.to_string().repeat(n)
}

/// The forward read of 50 and its mate, which overlap by 20.
fn forward() -> String {
    FRAGMENT[..50].to_string()
}

fn reverse() -> String {
    rc(&FRAGMENT[30..])
}

/// `seq` with the base at `i` replaced by another.
fn substitute(seq: &str, i: usize) -> String {
    let base = if &seq[i..=i] == "A" { "C" } else { "A" };
    format!("{}{base}{}", &seq[..i], &seq[i + 1..])
}

/// What a run wrote: its status and error, the merged pairs, the forward
/// and reverse reads not merged, the report and the table of errors.
struct Merged {
    output: Output,
    merged: String,
    not_merged: String,
    rev_not_merged: String,
    report: String,
    ee_table: String,
}

/// Runs mergpa on `forward`, fed on stdin, and `reverse`, from a file named
/// after `name`, with `args`, writing every output.
fn merge(name: &str, args: &[&str], forward: &str, reverse: &str) -> Merged {
    let paths = ["rev", "not-merged", "rev-not-merged", "report", "ee"]
        .map(|p| temp_path(&format!("{name}-{p}")));
    std::fs::write(&paths[0], reverse).unwrap();
    let path = |i: usize| paths[i].to_str().unwrap().to_string();
    let mut all = vec![
        "--reverse".to_string(),
        path(0),
        "--not-merged".to_string(),
        path(1),
        "--rev-not-merged".to_string(),
        path(2),
        "--report".to_string(),
        path(3),
        "--ee-table".to_string(),
        path(4),
    ];
    all.extend(args.iter().map(|a| a.to_string()));
    let mut child = Command::new(env!("CARGO_BIN_EXE_mergpa"))
        .args(&all)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("mergpa should start");
    child
        .stdin
        .take()
        .unwrap()
        .write_all(forward.as_bytes())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    // an output not written reads as empty
    let read = |i: usize| std::fs::read_to_string(&paths[i]).unwrap_or_default();
    let merged = Merged {
        merged: String::from_utf8(output.stdout.clone()).unwrap(),
        not_merged: read(1),
        rev_not_merged: read(2),
        report: read(3),
        ee_table: read(4),
        output,
    };
    for path in &paths {
        let _ = std::fs::remove_file(path);
    }
    merged
}

/// Runs mergpa, which must succeed.
fn merge_ok(name: &str, args: &[&str], forward: &str, reverse: &str) -> Merged {
    let merged = merge(name, args, forward, reverse);
    assert!(
        merged.output.status.success(),
        "mergpa failed: {}",
        String::from_utf8_lossy(&merged.output.stderr)
    );
    merged
}

/// Asserts that a run failed with exit code 1 and `message`.
fn assert_fails(merged: &Merged, message: &str) {
    assert_eq!(merged.output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&merged.output.stderr);
    assert!(stderr.contains(message), "{stderr}");
}

#[test]
fn merges_overlapping_pairs() {
    // a fragment of 50 read 30 from each end; bases are upper-cased, the
    // forward header kept whole, and agreeing bases scored at most Q41
    let fragment = "ACGTTGCAAGGCTTACCGATGCATCGGATCCTAGGTCAACGTTAGCCTAG";
    let (fwd, rev) = (&fragment[..30], rc(&fragment[20..]));
    let forward = fastq(&[
        ("a/1 x", fwd, &q('I', 30)),
        ("b", &fwd.to_lowercase(), &q('5', 30)),
    ]);
    let reverse = fastq(&[("a/2", &rev, &q('5', 30)), ("b2", &rev, &q('I', 30))]);
    let merged = merge_ok("basic", &[], &forward, &reverse);
    assert_eq!(
        merged.merged,
        format!(
            "@a/1 x\n{fragment}\n+\n{}{}{}\n@b\n{fragment}\n+\n{}{}{}\n",
            q('I', 20),
            q('J', 10),
            q('5', 20),
            q('5', 20),
            q('J', 10),
            q('I', 20),
        )
    );
    assert_eq!(merged.not_merged, "");
    assert_eq!(merged.rev_not_merged, "");
    // expected errors of the overlap in each read, and differences
    assert_eq!(
        merged.ee_table,
        "0.003000\t0.3000\t0\t0\n0.3000\t0.003000\t0\t0\n"
    );
    assert_eq!(
        merged.report,
        "         2  Pairs
         2  Merged (100.0%)
         0  Not merged (0.0%)

Statistics of all reads:
     30.00  Mean read length

Statistics of merged reads:
     50.00  Mean fragment length
      0.00  Standard deviation of fragment length
      0.15  Mean expected error in forward sequences
      0.15  Mean expected error in reverse sequences
      0.20  Mean expected error in merged sequences
      0.00  Mean observed errors in merged region of forward sequences
      0.00  Mean observed errors in merged region of reverse sequences
      0.00  Mean observed errors in merged region
"
    );
}

#[test]
fn resolves_differences_and_ns() {
    // a difference goes to the higher score, with the posterior score of a
    // disagreement; an N takes the other read's base and score
    let changed = substitute(&forward(), 35);
    let forward_n = format!("{}N{}", &FRAGMENT[..40], &FRAGMENT[41..50]);
    let reverse_n = format!("{}N{}", &reverse()[..5], &reverse()[6..]);
    let fwd = fastq(&[
        ("hi", &changed, &format!("{}5{}", q('I', 35), q('I', 14))),
        ("lo", &changed, &q('5', 50)),
        ("nf", &forward_n, &q('I', 50)),
        ("nr", &forward(), &q('I', 50)),
    ]);
    let rev = fastq(&[
        ("hi", &reverse(), &q('5', 50)),
        ("lo", &reverse(), &q('I', 50)),
        ("nf", &reverse(), &q('I', 50)),
        ("nr", &reverse_n, &q('I', 50)),
    ]);
    let merged = merge_ok("diffs", &[], &fwd, &rev);
    let with_n = format!("{}N{}", &FRAGMENT[..74], &FRAGMENT[75..]);
    assert_eq!(
        merged.merged,
        format!(
            "@hi\n{FRAGMENT}\n+\n{}JJJJJ$JJJJJJJJJJJJJJ{}\n\
             @lo\n{FRAGMENT}\n+\n{}JJJJJ5JJJJJJJJJJJJJJ{}\n\
             @nf\n{FRAGMENT}\n+\n{}JJJJJJJJJJIJJJJJJJJJ{}\n\
             @nr\n{with_n}\n+\n{}{}{}!IIIII\n",
            q('I', 30),
            q('5', 30),
            q('5', 30),
            q('I', 30),
            q('I', 30),
            q('I', 30),
            q('I', 30),
            q('J', 20),
            q('I', 24),
        )
    );
    assert_eq!(
        merged.ee_table,
        "0.01490\t0.5000\t1\t0\n0.5000\t0.005000\t1\t0\n\
         0.7549\t0.005000\t1\t0\n0.005000\t0.7549\t0\t0\n"
    );
    assert!(
        merged
            .report
            .contains("      0.75  Mean observed errors in merged region of forward sequences\n")
    );
}

#[test]
fn bounds_merged_scores() {
    let fwd = fastq(&[("a", &forward(), &q('I', 50))]);
    let rev = fastq(&[("a", &reverse(), &q('I', 50))]);
    let merged = merge_ok("qout", &["--qmax-out", "93", "--qmin-out", "5"], &fwd, &rev);
    assert_eq!(
        merged.merged,
        format!(
            "@a\n{FRAGMENT}\n+\n{}{}{}\n",
            q('I', 30),
            q('v', 20),
            q('I', 30)
        )
    );
}

#[test]
fn reads_scores_at_offset_64() {
    let fwd = fastq(&[("a", &forward(), &q('h', 50))]);
    let rev = fastq(&[("a", &reverse(), &q('T', 50))]);
    let merged = merge_ok("ascii", &["--ascii", "64"], &fwd, &rev);
    assert_eq!(
        merged.merged,
        format!(
            "@a\n{FRAGMENT}\n+\n{}{}{}\n",
            q('h', 30),
            q('i', 20),
            q('T', 30)
        )
    );
    assert_eq!(merged.ee_table, "0.005000\t0.5000\t0\t0\n");
}

#[test]
fn merges_staggered_pairs_with_stagger() {
    // the reverse read runs past the start of the forward read
    let fwd = fastq(&[("s", &FRAGMENT[10..60], &q('I', 50))]);
    let rev = fastq(&[("s", &rc(&FRAGMENT[5..55]), &q('5', 50))]);
    let merged = merge_ok("staggered", &[], &fwd, &rev);
    assert_eq!(merged.merged, "");
    assert_eq!(merged.not_merged, fwd);
    assert_eq!(merged.rev_not_merged, rev);
    assert!(merged.report.contains("         1  staggered read pairs\n"));
    // merged, the overhangs are cut
    let merged = merge_ok("stagger", &["--stagger"], &fwd, &rev);
    assert_eq!(
        merged.merged,
        format!("@s\n{}\n+\n{}\n", &FRAGMENT[10..55], q('J', 45))
    );
    assert_eq!(merged.ee_table, "0.004500\t0.4500\t0\t0\n");
}

#[test]
fn limits_differences() {
    // three differences in an overlap of 36
    let fragment = long_fragment();
    let changed = substitute(&substitute(&substitute(&fragment[..70], 40), 50), 60);
    let fwd = fastq(&[("t", &changed, &q('I', 70))]);
    let rev = fastq(&[("t", &rc(&fragment[34..]), &q('I', 70))]);
    let expected = format!(
        "@t\n{fragment}\n+\n{}JJJJJJ$JJJJJJJJJ$JJJJJJJJJ$JJJJJJJJJ{}\n",
        q('I', 34),
        q('I', 34)
    );
    for (args, merges, reason) in [
        (&["--max-diffs", "2"][..], false, "too many differences"),
        (&["--max-diffs", "3"], true, ""),
        (
            &["--max-diff-pct", "8"],
            false,
            "too high percentage of differences",
        ),
        (&["--max-diff-pct", "9"], true, ""),
    ] {
        let merged = merge_ok("diffs", args, &fwd, &rev);
        if merges {
            assert_eq!(merged.merged, expected, "{args:?}");
            assert_eq!(merged.ee_table, "0.007000\t0.007000\t3\t0\n");
        } else {
            assert_eq!(merged.merged, "", "{args:?}");
            assert_eq!(merged.not_merged, fwd);
            assert!(merged.report.contains(reason), "{args:?}");
        }
    }
}

#[test]
fn requires_the_minimum_overlap() {
    // an overlap of 9 is too short by default
    let fwd = fastq(&[("o", &FRAGMENT[..40], &q('I', 40))]);
    let rev = fastq(&[("o", &rc(&FRAGMENT[31..71]), &q('I', 40))]);
    let merged = merge_ok("ovl9", &[], &fwd, &rev);
    assert_eq!(merged.merged, "");
    assert!(merged.report.contains("         1  overlap too short\n"));
    let merged = merge_ok("ovl9b", &["--min-overlap", "9"], &fwd, &rev);
    assert_eq!(
        merged.merged,
        format!(
            "@o\n{}\n+\n{}{}{}\n",
            &FRAGMENT[..71],
            q('I', 31),
            q('J', 9),
            q('I', 31)
        )
    );
    // below 9 the k-mer and score thresholds relax
    let fwd = fastq(&[("o", &FRAGMENT[..40], &q('I', 40))]);
    let rev = fastq(&[("o", &rc(&FRAGMENT[33..73]), &q('I', 40))]);
    assert_eq!(merge_ok("ovl7", &[], &fwd, &rev).merged, "");
    assert_eq!(
        merge_ok("ovl7-8", &["--min-overlap", "8"], &fwd, &rev).merged,
        ""
    );
    let merged = merge_ok("ovl7-5", &["--min-overlap", "5"], &fwd, &rev);
    assert_eq!(
        merged.merged,
        format!(
            "@o\n{}\n+\n{}{}{}\n",
            &FRAGMENT[..73],
            q('I', 33),
            q('J', 7),
            q('I', 33)
        )
    );
}

#[test]
fn reports_pairs_not_merged() {
    // in order: two merged pairs of different lengths, an empty read, a
    // repeat, low scores, and a staggered pair. Reads not merged are
    // written upper-cased, an N scored 0
    let repeat = "ACGTGCAT".repeat(10);
    let low = format!("{}N{}", FRAGMENT[..20].to_lowercase(), &FRAGMENT[21..50]);
    let fwd = fastq(&[
        ("m1", &forward(), &q('I', 50)),
        ("m2", &FRAGMENT[..40], &q('I', 40)),
        ("e", "", ""),
        ("rep", &repeat[..50], &q('I', 50)),
        ("low", &low, &q('"', 50)),
        ("st", &FRAGMENT[10..60], &q('I', 50)),
    ]);
    let rev = fastq(&[
        ("m1", &reverse(), &q('I', 50)),
        ("m2", &rc(&FRAGMENT[25..70]), &q('I', 45)),
        ("e", &reverse(), &q('I', 50)),
        ("rep", &rc(&repeat[30..80]), &q('I', 50)),
        ("low", &reverse(), &q('#', 50)),
        ("st", &rc(&FRAGMENT[5..55]), &q('I', 50)),
    ]);
    let merged = merge_ok("mixed", &[], &fwd, &rev);
    assert_eq!(
        merged.merged,
        format!(
            "@m1\n{FRAGMENT}\n+\n{}{}{}\n@m2\n{}\n+\n{}{}{}\n",
            q('I', 30),
            q('J', 20),
            q('I', 30),
            &FRAGMENT[..70],
            q('I', 25),
            q('J', 15),
            q('I', 30)
        )
    );
    assert_eq!(
        merged.not_merged,
        format!(
            "@e\n\n+\n\n@rep\n{}\n+\n{}\n@low\n{}N{}\n+\n{}!{}\n@st\n{}\n+\n{}\n",
            &repeat[..50],
            q('I', 50),
            &FRAGMENT[..20],
            &FRAGMENT[21..50],
            q('"', 20),
            q('"', 29),
            &FRAGMENT[10..60],
            q('I', 50)
        )
    );
    assert_eq!(
        merged.rev_not_merged,
        fastq(&[
            ("e", &reverse(), &q('I', 50)),
            ("rep", &rc(&repeat[30..80]), &q('I', 50)),
            ("low", &reverse(), &q('#', 50)),
            ("st", &rc(&FRAGMENT[5..55]), &q('I', 50)),
        ])
    );
    assert_eq!(
        merged.ee_table,
        "0.005000\t0.005000\t0\t0\n0.004000\t0.004500\t0\t0\n"
    );
    assert_eq!(
        merged.report,
        "         6  Pairs
         2  Merged (33.3%)
         4  Not merged (66.7%)

Pairs that failed merging due to various reasons:
         1  reads too short (after truncation)
         1  multiple potential alignments
         1  alignment score too low, or score drop too high
         1  staggered read pairs

Statistics of all reads:
     44.58  Mean read length

Statistics of merged reads:
     75.00  Mean fragment length
      5.00  Standard deviation of fragment length
      0.00  Mean expected error in forward sequences
      0.00  Mean expected error in reverse sequences
      0.01  Mean expected error in merged sequences
      0.00  Mean observed errors in merged region of forward sequences
      0.00  Mean observed errors in merged region of reverse sequences
      0.00  Mean observed errors in merged region
"
    );
}

#[test]
fn reports_pairs_without_shared_kmers() {
    let reads = fastq(&[("a", "ACGT", "IIII")]);
    let merged = merge_ok("nokmers", &[], &reads, &reads);
    assert_eq!(merged.not_merged, reads);
    assert_eq!(merged.rev_not_merged, reads);
    assert_eq!(
        merged.report,
        "         1  Pairs
         0  Merged (0.0%)
         1  Not merged (100.0%)

Pairs that failed merging due to various reasons:
         1  too few kmers found on same diagonal

Statistics of all reads:
      4.00  Mean read length
"
    );
}

#[test]
fn merges_empty_input() {
    let merged = merge_ok("nothing", &[], "", "");
    assert_eq!(merged.merged, "");
    assert_eq!(merged.ee_table, "");
    assert_eq!(
        merged.report,
        "         0  Pairs\n         0  Merged\n         0  Not merged\n\n"
    );
}

#[test]
fn checks_the_range_of_scores() {
    // the scores of a pair with an empty read are not checked
    let fwd = fastq(&[("e", "", ""), ("f", &forward(), &q('h', 50))]);
    let rev = fastq(&[("e", &reverse(), &q('h', 50)), ("f", "", "")]);
    let merged = merge_ok("unchecked", &["--qmax", "41"], &fwd, &rev);
    assert_eq!(merged.not_merged, fwd);
    // h is 71
    let fwd = fastq(&[
        ("e", &forward(), &q('I', 50)),
        ("f", &forward(), &q('h', 50)),
    ]);
    let rev = fastq(&[
        ("e", &reverse(), &q('I', 50)),
        ("f", &reverse(), &q('I', 50)),
    ]);
    assert_fails(&merge("above", &["--qmax", "41"], &fwd, &rev), "above");
    let fwd = fastq(&[
        ("e", &forward(), &q('I', 50)),
        ("f", &forward(), &q('I', 50)),
    ]);
    let rev = fastq(&[
        ("e", &reverse(), &q('I', 50)),
        ("f", &reverse(), &format!("{}!", q('I', 49))),
    ]);
    assert_fails(&merge("below", &["--qmin", "1"], &fwd, &rev), "below");
}

#[test]
fn checks_abundances_of_both_reads() {
    let reads = fastq(&[("a", &forward(), &q('I', 50))]);
    let zero = fastq(&[("a;size=0", &reverse(), &q('I', 50))]);
    assert_fails(&merge("rev-zero", &[], &reads, &zero), "abundance");
    let zero = fastq(&[("a;size=0", &forward(), &q('I', 50))]);
    let reads = fastq(&[("a", &reverse(), &q('I', 50))]);
    assert_fails(&merge("fwd-zero", &[], &zero, &reads), "abundance");
}

#[test]
fn fails_on_unequal_numbers_of_reads() {
    let one = fastq(&[("a", &forward(), &q('I', 50))]);
    let two = fastq(&[
        ("a", &forward(), &q('I', 50)),
        ("b", &forward(), &q('I', 50)),
    ]);
    assert_fails(
        &merge("more-forward", &[], &two, &one),
        "More forward reads than reverse reads",
    );
    assert_fails(
        &merge("more-reverse", &[], &one, &two),
        "More reverse reads than forward reads",
    );
    assert_fails(
        &merge("no-forward", &[], "", &one),
        "More reverse reads than forward reads",
    );
}

#[test]
fn rejects_fasta() {
    let reads = fastq(&[("a", &forward(), &q('I', 50))]);
    let fasta = ">a\nAC\n";
    for (forward, reverse) in [
        (fasta, reads.as_str()),
        (reads.as_str(), fasta),
        ("", fasta),
    ] {
        assert_fails(
            &merge("fasta", &[], forward, reverse),
            "FASTQ file expected",
        );
    }
}

#[test]
fn rejects_invalid_input() {
    let reads = fastq(&[("a", &forward(), &q('I', 50))]);
    let merged = merge("invalid", &[], "@a\nACGT\nIIII\n", &reads);
    assert_eq!(merged.output.status.code(), Some(1));
    assert!(merged.merged.is_empty());
    let merged = merge("invalid-rev", &[], &reads, "@a\nAC\n+\nI\n");
    assert_eq!(merged.output.status.code(), Some(1));
}

#[test]
fn rejects_invalid_arguments() {
    // vsearch: an overlap of at least 5, differences of at least 0, a
    // percentage from 0 to 100, and score bounds within the printable range
    let reads = fastq(&[("a", &forward(), &q('I', 50))]);
    for args in [
        &["--unknown"][..],
        &["--min-overlap", "4"],
        &["--max-diffs", "-1"],
        &["--max-diff-pct", "-1"],
        &["--max-diff-pct", "100.5"],
        &["--ascii", "50"],
        &["--qmin", "10", "--qmax", "5"],
        &["--ascii", "64", "--qmax", "63"],
        &["--qmin", "-1"],
        &["--qmin-out", "42"],
        &["--qmin-out", "10", "--qmax-out", "5"],
        &["--qmax-out", "94"],
        &["--ascii", "64", "--qmax-out", "63"],
        &["--qmin-out", "-1"],
    ] {
        let merged = merge("arguments", args, &reads, &reads);
        assert_eq!(merged.output.status.code(), Some(2), "{args:?} should fail");
        assert!(merged.merged.is_empty());
    }
    // --reverse is required
    let output = Command::new(env!("CARGO_BIN_EXE_mergpa"))
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn reads_and_writes_files() {
    let (input, rev, out) = (temp_path("in.fq"), temp_path("rev.fq"), temp_path("out.fq"));
    std::fs::write(&input, fastq(&[("a", &forward(), &q('I', 50))])).unwrap();
    std::fs::write(&rev, fastq(&[("a", &reverse(), &q('I', 50))])).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_mergpa"))
        .args([
            "-i",
            input.to_str().unwrap(),
            "--reverse",
            rev.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    for path in [input, rev, out] {
        std::fs::remove_file(path).unwrap();
    }
    assert!(output.status.success());
    // nothing on stdout or stderr without --report
    assert!(output.stdout.is_empty());
    assert!(output.stderr.is_empty());
    assert!(
        written.starts_with(&format!("@a\n{FRAGMENT}\n")),
        "{written}"
    );
}
