//! End-to-end tests of the syncpa binary.
//!
//! Expected outputs are those of `vsearch --fastx_syncpairs --fasta_width
//! 0` with `--reverse` on the same inputs. Forward reads whose mate is in
//! the reverse file are written in forward order to `-o` (`--fastaout`,
//! `--fastqout`), and their mates in the same order to `--rev-out`
//! (`--fastaout_rev`). Forward reads without a mate go to `--orphans`
//! (`--fastaout_orphans`). Reverse reads without one go, in reverse-file
//! order, to `--rev-orphans` (`--fastaout_orphans_rev`). Mates share their
//! header up to the first space or tab, less a final `/1` or `/2`. The
//! characters that may come before the 1 or 2 are `--separators`
//! (`--read_separators`, default `/`).

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-syncpa-{}-{name}", std::process::id()))
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

/// What a run wrote: its status and error, and the synced forward and
/// reverse reads and the forward and reverse orphans.
struct Synced {
    output: Output,
    forward: String,
    reverse: String,
    orphans: String,
    rev_orphans: String,
}

/// Runs syncpa on `forward`, fed on stdin, and `reverse`, from a file named
/// after `name`, with `args`, writing every output.
fn sync(name: &str, args: &[&str], forward: &str, reverse: &str) -> Synced {
    let paths =
        ["rev", "rev-out", "orphans", "rev-orphans"].map(|p| temp_path(&format!("{name}-{p}")));
    std::fs::write(&paths[0], reverse).unwrap();
    let path = |i: usize| paths[i].to_str().unwrap().to_string();
    let mut all = vec![
        "--reverse".to_string(),
        path(0),
        "--rev-out".to_string(),
        path(1),
        "--orphans".to_string(),
        path(2),
        "--rev-orphans".to_string(),
        path(3),
    ];
    all.extend(args.iter().map(|a| a.to_string()));
    let mut child = Command::new(env!("CARGO_BIN_EXE_syncpa"))
        .args(&all)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("syncpa should start");
    let mut stdin = child.stdin.take().unwrap();
    let fed = forward.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    // an output not written reads as empty
    let read = |i: usize| std::fs::read_to_string(&paths[i]).unwrap_or_default();
    let synced = Synced {
        forward: String::from_utf8(output.stdout.clone()).unwrap(),
        reverse: read(1),
        orphans: read(2),
        rev_orphans: read(3),
        output,
    };
    for path in &paths {
        let _ = std::fs::remove_file(path);
    }
    synced
}

/// Runs syncpa as `sync` does, which must succeed.
fn synced(name: &str, args: &[&str], forward: &str, reverse: &str) -> Synced {
    let synced = sync(name, args, forward, reverse);
    assert!(
        synced.output.status.success(),
        "syncpa failed: {}",
        String::from_utf8_lossy(&synced.output.stderr)
    );
    synced
}

/// Forward reads with mates by /1 and /2, by Casava descriptions, by _1
/// and _2, and none; a header with an annotation of zero.
fn forward() -> String {
    fastq(&[
        ("a/1", "ACGT", "IIII"),
        ("b 1:N", "CC", "II"),
        ("c/1", "GG", "#I"),
        ("d_1", "TT", "II"),
        ("e;size=0/1", "AA", "II"),
    ])
}

/// The reverse reads of `forward`, in another order.
fn reverse() -> String {
    fastq(&[
        ("e;size=0/2", "TT", "II"),
        ("x/2", "GA", "II"),
        ("b 2:N", "GG", "5I"),
        ("a/2", "TTTT", "III#"),
        ("d_2", "AA", "II"),
    ])
}

#[test]
fn pairs_mates_by_label() {
    // abundances are not read, and a zero is no error
    let out = synced("pairs", &[], &forward(), &reverse());
    assert_eq!(
        out.forward,
        fastq(&[
            ("a/1", "ACGT", "IIII"),
            ("b 1:N", "CC", "II"),
            ("e;size=0/1", "AA", "II"),
        ])
    );
    assert_eq!(
        out.reverse,
        fastq(&[
            ("a/2", "TTTT", "III#"),
            ("b 2:N", "GG", "5I"),
            ("e;size=0/2", "TT", "II"),
        ])
    );
    assert_eq!(
        out.orphans,
        fastq(&[("c/1", "GG", "#I"), ("d_1", "TT", "II")])
    );
    assert_eq!(
        out.rev_orphans,
        fastq(&[("x/2", "GA", "II"), ("d_2", "AA", "II")])
    );
}

#[test]
fn strips_the_separators_given() {
    let out = synced(
        "separators",
        &["--separators", "_/"],
        &forward(),
        &reverse(),
    );
    assert_eq!(
        out.forward,
        fastq(&[
            ("a/1", "ACGT", "IIII"),
            ("b 1:N", "CC", "II"),
            ("d_1", "TT", "II"),
            ("e;size=0/1", "AA", "II"),
        ])
    );
    assert_eq!(
        out.reverse,
        fastq(&[
            ("a/2", "TTTT", "III#"),
            ("b 2:N", "GG", "5I"),
            ("d_2", "AA", "II"),
            ("e;size=0/2", "TT", "II"),
        ])
    );
    assert_eq!(out.orphans, fastq(&[("c/1", "GG", "#I")]));
    assert_eq!(out.rev_orphans, fastq(&[("x/2", "GA", "II")]));
}

#[test]
fn matches_up_to_the_first_blank() {
    // a key may be empty; /3 is kept; a tab ends the key; a single
    // character is never stripped; case is kept and FASTA is unwrapped
    let long = "A".repeat(100);
    let forward = fasta(&[
        ("/1", "AC"),
        ("a/3", "AC"),
        ("p\tq/1", "ac"),
        ("x/1 y", "GT"),
        ("1", "A"),
        ("2", &long),
    ]);
    let reverse = fasta(&[
        ("/2", "GG"),
        ("a/3", "TT"),
        ("p\tz/2", "CC"),
        ("x/2 z", "AA"),
        ("1", "T"),
        ("2/2", "G"),
    ]);
    let out = synced("blank", &[], &forward, &reverse);
    assert_eq!(out.forward, forward);
    assert_eq!(out.reverse, reverse);
    assert_eq!(out.orphans, "");
    assert_eq!(out.rev_orphans, "");
}

#[test]
fn fails_on_duplicate_reverse_labels() {
    // "/2" and "" share the empty key; nothing is written
    let out = sync(
        "dup-rev",
        &[],
        &fasta(&[("a/1", "A")]),
        &fasta(&[("/2", "GG"), ("", "G")]),
    );
    assert_eq!(out.output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&out.output.stderr);
    assert!(
        stderr.contains("Duplicate read label in reverse file"),
        "{stderr}"
    );
    assert_eq!(
        out.forward + &out.reverse + &out.orphans + &out.rev_orphans,
        ""
    );
}

#[test]
fn fails_on_duplicate_forward_labels_that_pair() {
    // the reads before the second mate are written
    let out = sync(
        "dup-fwd",
        &[],
        &fasta(&[("a/1", "A"), ("b/1", "C"), ("a/1", "G")]),
        &fasta(&[("a/2", "T")]),
    );
    assert_eq!(out.output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&out.output.stderr);
    assert!(
        stderr.contains("Duplicate read label in forward file"),
        "{stderr}"
    );
    assert_eq!(out.forward, fasta(&[("a/1", "A")]));
    assert_eq!(out.reverse, fasta(&[("a/2", "T")]));
    assert_eq!(out.orphans, fasta(&[("b/1", "C")]));
    assert_eq!(out.rev_orphans, "");
}

#[test]
fn writes_forward_orphans_sharing_a_label() {
    let out = synced(
        "dup-orphans",
        &[],
        &fasta(&[("b/1", "A"), ("b/1", "C")]),
        &fasta(&[("a/2", "T")]),
    );
    assert_eq!(out.forward, "");
    assert_eq!(out.orphans, fasta(&[("b/1", "A"), ("b/1", "C")]));
    assert_eq!(out.rev_orphans, fasta(&[("a/2", "T")]));
}

#[test]
fn takes_the_format_of_the_file_with_reads() {
    let out = synced("empty-fwd", &[], "", &reverse());
    assert_eq!(out.forward + &out.reverse + &out.orphans, "");
    assert_eq!(out.rev_orphans, reverse());
    let out = synced("empty-rev", &[], &forward(), "");
    assert_eq!(out.orphans, forward());
    let out = synced("empty", &[], "", "");
    assert_eq!(
        out.forward + &out.reverse + &out.orphans + &out.rev_orphans,
        ""
    );
}

#[test]
fn rejects_mixed_formats() {
    let out = sync("mixed", &[], &fasta(&[("a/1", "A")]), &reverse());
    assert_eq!(out.output.status.code(), Some(1));
    assert_eq!(
        out.forward + &out.reverse + &out.orphans + &out.rev_orphans,
        ""
    );
}

#[test]
fn rejects_invalid_input() {
    // in the reverse file, before anything is written
    let out = sync(
        "bad-rev",
        &[],
        &fasta(&[("a/1", "A")]),
        &fasta(&[("a/2", "T"), ("b/2", "A-C")]),
    );
    assert_eq!(out.output.status.code(), Some(1));
    assert_eq!(
        out.forward + &out.reverse + &out.orphans + &out.rev_orphans,
        ""
    );
    // in the forward file, after the pairs before it
    let out = sync(
        "bad-fwd",
        &[],
        &fasta(&[("a/1", "A"), ("b/1", "A-C")]),
        &fasta(&[("a/2", "T")]),
    );
    assert_eq!(out.output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&out.output.stderr);
    assert!(
        stderr.contains("Illegal sequence character '-'"),
        "{stderr}"
    );
    assert_eq!(out.forward, fasta(&[("a/1", "A")]));
    assert_eq!(out.reverse, fasta(&[("a/2", "T")]));
    assert_eq!(out.orphans + &out.rev_orphans, "");
}

#[test]
fn rejects_invalid_arguments() {
    for args in [&[][..], &["--unknown"], &["--rev-out", "x"]] {
        let output = Command::new(env!("CARGO_BIN_EXE_syncpa"))
            .args(args)
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, rev, output) = (temp_path("in.fq"), temp_path("r.fq"), temp_path("out.fq"));
    std::fs::write(&input, forward()).unwrap();
    std::fs::write(&rev, reverse()).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_syncpa"))
        .args([
            "-i",
            input.to_str().unwrap(),
            "--reverse",
            rev.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&output).unwrap();
    for path in [input, rev, output] {
        std::fs::remove_file(path).unwrap();
    }
    assert!(status.status.success());
    assert!(status.stdout.is_empty());
    assert_eq!(
        written,
        fastq(&[
            ("a/1", "ACGT", "IIII"),
            ("b 1:N", "CC", "II"),
            ("e;size=0/1", "AA", "II"),
        ])
    );
}
