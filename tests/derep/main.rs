//! End-to-end tests of the derep binary.
//!
//! Unless noted otherwise, expected outputs are those of `vsearch
//! --fastx_uniques --sizein --sizeout --minseqlength 0 --notrunclabels
//! --fasta_width 0` on the same input. derep keeps whole headers and every
//! length, and leaves filtering, relabelling and format changes to other
//! tools.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs derep with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_derep"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("derep should start");
    let mut stdin = child.stdin.take().unwrap();
    let fed = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// Runs derep and returns its standard output, which must succeed.
fn derep(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "derep failed: {}",
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
    std::env::temp_dir().join(format!("rsearch-derep-{}-{name}", std::process::id()))
}

#[test]
fn merges_identical_sequences() {
    let input = fasta(&[("a", "ACGT"), ("b", "ACGT"), ("c", "TTTT"), ("d", "acgt")]);
    let expected = fasta(&[("a;size=3", "ACGT"), ("c;size=1", "TTTT")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn ignores_case_and_u_versus_t() {
    // the first sequence is written as it appears in the input
    let input = fasta(&[("a", "ACGU"), ("b", "acgt")]);
    assert_eq!(derep(&[], &input), fasta(&[("a;size=2", "ACGU")]));
}

#[test]
fn keeps_ambiguity_codes_distinct() {
    let input = fasta(&[("a", "ACGTN"), ("b", "ACGTR"), ("c", "ACGTn")]);
    let expected = fasta(&[("a;size=2", "ACGTN"), ("b;size=1", "ACGTR")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn adds_up_input_sizes() {
    let input = fasta(&[("a;size=2", "AAAA"), ("b;size=3", "CCCC"), ("c", "AAAA")]);
    let expected = fasta(&[("a;size=3", "AAAA"), ("b;size=3", "CCCC")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn sorts_by_decreasing_size_then_label() {
    let input = fasta(&[
        ("b", "CCCC"),
        ("a", "GGGG"),
        ("c", "AAAA"),
        ("d", "TTTT"),
        ("e", "TTTT"),
    ]);
    let expected = fasta(&[
        ("d;size=2", "TTTT"),
        ("a;size=1", "GGGG"),
        ("b;size=1", "CCCC"),
        ("c;size=1", "AAAA"),
    ]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn compares_whole_labels_to_break_ties() {
    let input = fasta(&[("b;x=1", "AAAA"), ("b", "CCCC"), ("b;a=1", "GGGG")]);
    let expected = fasta(&[
        ("b;size=1", "CCCC"),
        ("b;a=1;size=1", "GGGG"),
        ("b;x=1;size=1", "AAAA"),
    ]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn breaks_ties_on_label_by_first_appearance() {
    let input = fasta(&[("x", "CCCC"), ("x", "AAAA")]);
    let expected = fasta(&[("x;size=1", "CCCC"), ("x;size=1", "AAAA")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn unsorted_keeps_order_of_first_appearance() {
    // vsearch --derep_smallmem
    let input = fasta(&[("d", "TTTT"), ("c", "GGGG"), ("a", "TTTT")]);
    let expected = fasta(&[("d;size=2", "TTTT"), ("c;size=1", "GGGG")]);
    assert_eq!(derep(&["--unsorted"], &input), expected);
}

#[test]
fn merges_reverse_complements_only_with_strand() {
    let input = fasta(&[("a", "AACC"), ("b", "GGTT"), ("c", "GGTT")]);
    let plus = fasta(&[("b;size=2", "GGTT"), ("a;size=1", "AACC")]);
    let both = fasta(&[("a;size=3", "AACC")]);
    assert_eq!(derep(&[], &input), plus);
    assert_eq!(derep(&["--strand"], &input), both);
    assert_eq!(derep(&["-s"], &input), both);
}

#[test]
fn prefers_the_plus_strand_match() {
    // GGTT is both a copy of the second sequence and the reverse
    // complement of the first; the identical sequence wins
    let input = fasta(&[("a", "AACC"), ("b", "GGTT"), ("c", "GGTT")]);
    let output = derep(&["--unsorted", "--strand"], &input);
    assert_eq!(output, fasta(&[("a;size=3", "AACC")]));
    let input = fasta(&[("b", "GGTT"), ("a", "AACC"), ("c", "GGTT")]);
    let output = derep(&["--unsorted", "--strand"], &input);
    assert_eq!(output, fasta(&[("b;size=3", "GGTT")]));
}

#[test]
fn keeps_empty_sequences() {
    // vsearch --minseqlength 0: empty sequences are dereplicated too
    let input = fasta(&[("a", "AAAA"), ("e", ""), ("b;size=3", "CCCC"), ("f", "")]);
    let expected = fasta(&[("b;size=3", "CCCC"), ("e;size=2", ""), ("a;size=1", "AAAA")]);
    assert_eq!(derep(&[], &input), expected);
    assert_eq!(derep(&["--strand"], &input), expected);
    // vsearch --derep_smallmem
    let expected = fasta(&[("a;size=1", "AAAA"), ("e;size=2", ""), ("b;size=3", "CCCC")]);
    assert_eq!(derep(&["--unsorted"], &input), expected);
}

#[test]
fn joins_multiline_sequences() {
    let input = ">a\nAC\nGT\n>b\nACGT\n";
    assert_eq!(derep(&[], input), fasta(&[("a;size=2", "ACGT")]));
}

#[test]
fn replaces_only_the_size_annotation() {
    // other annotations and descriptions are kept, a size followed by a
    // description is not read, and the first valid size annotation is
    // replaced by the new size at the end
    for (label, expected) in [
        (
            "sample=1;foo;size=3 desc",
            "sample=1;foo;size=3 desc;size=2",
        ),
        ("a;size=3;x=1", "a;x=1;size=4"),
        ("a;size=2;size=5", "a;size=5;size=3"),
        ("size=3", ";size=4"),
        ("a;size=abc", "a;size=abc;size=2"),
        ("asize=2", "asize=2;size=2"),
        ("a;b;size=2 x", "a;b;size=2 x;size=2"),
    ] {
        let input = fasta(&[(label, "ACGT"), ("z", "ACGT")]);
        assert_eq!(derep(&[], &input), fasta(&[(expected, "ACGT")]), "{label}");
    }
}

#[test]
fn rejects_zero_sizes() {
    let output = run(&[], &fasta(&[("a;size=0", "ACGT")]));
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("zero"));
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fa"), temp_path("out.fa"));
    std::fs::write(&input, fasta(&[("a", "ACGT"), ("b", "ACGT")])).unwrap();
    let stdout = derep(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ],
        "",
    );
    assert_eq!(stdout, "");
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(written, fasta(&[("a;size=2", "ACGT")]));
}

#[test]
fn merges_fastq_qualities() {
    let input = fastq(&[("a", "ACGT", "I#5+"), ("b", "ACGT", "#I5,")]);
    assert_eq!(derep(&[], &input), fastq(&[("a;size=2", "ACGT", "&&5+")]));
}

#[test]
fn weights_fastq_qualities_by_size() {
    let input = fastq(&[
        ("a;size=3", "ACGT", "5555"),
        ("b;size=2", "acgu", "++++"),
        ("c", "ACGT", "????"),
    ]);
    assert_eq!(derep(&[], &input), fastq(&[("a;size=6", "ACGT", "....")]));
}

#[test]
fn reverses_qualities_of_reverse_strand_reads() {
    // vsearch --fastx_uniques --strand both
    let input = fastq(&[("a", "AACC", "ABCD"), ("b", "GGTT", "IIII")]);
    let expected = fastq(&[("a;size=2", "AACC", "CDEE")]);
    assert_eq!(derep(&["--strand"], &input), expected);
}

#[test]
fn keeps_empty_fastq_reads() {
    let input = fastq(&[("a", "ACGT", "IIII"), ("e", "", ""), ("b", "ACGT", "IIII")]);
    let expected = fastq(&[("a;size=2", "ACGT", "IIII"), ("e;size=1", "", "")]);
    assert_eq!(derep(&[], &input), expected);
}

#[test]
fn keeps_whole_headers() {
    // vsearch --notrunclabels: a header goes on after a space or a tab, and
    // ends at a carriage return
    let input = ">a\tdesc\nACGT\n>b desc\nACGT\n";
    assert_eq!(derep(&[], input), fasta(&[("a\tdesc;size=2", "ACGT")]));
    let input = "@a\tdesc\nACGT\n+\nIIII\n@b\rc\nACGT\n+\nIIII\n";
    assert_eq!(
        derep(&[], input),
        fastq(&[("a\tdesc;size=2", "ACGT", "IIII")])
    );
}

#[test]
fn strips_invalid_fasta_characters_with_a_warning() {
    let output = run(&[], ">a\nAC1T\n>b\nACT\n");
    assert!(output.status.success());
    assert_eq!(output.stdout, fasta(&[("a;size=2", "ACT")]).as_bytes());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("1 invalid characters stripped from FASTA file: 1(1)"));
}

#[test]
fn fails_on_invalid_input() {
    for (input, message) in [
        ("xyz", "File type not recognized"),
        ("@a\nAC*T\n+\nIIII\n", "Illegal sequence character"),
        (">a\nAC-T\n", "Illegal sequence character"),
    ] {
        let output = run(&[], input);
        assert_eq!(output.status.code(), Some(1), "{input:?} should fail");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{input:?}: {stderr}");
    }
}

#[test]
fn raises_quality_zero_to_one() {
    let input = fastq(&[("a", "ACGT", "!!+I")]);
    assert_eq!(derep(&[], &input), fastq(&[("a;size=1", "ACGT", "\"\"+I")]));
}

#[test]
fn counts_each_read_once_without_sizein() {
    // vsearch without --sizein
    let input = fasta(&[("a;size=3", "ACGT"), ("b;size=2", "ACGT"), ("c", "TTTT")]);
    let expected = fasta(&[("a;size=2", "ACGT"), ("c;size=1", "TTTT")]);
    assert_eq!(derep(&["--no-sizein"], &input), expected);
    // a zero size is only refused when sizes are read
    let input = fasta(&[("a;size=0", "ACGT")]);
    assert_eq!(
        derep(&["--no-sizein"], &input),
        fasta(&[("a;size=1", "ACGT")])
    );
}

#[test]
fn keeps_labels_without_sizeout() {
    // vsearch without --sizeout
    let input = fasta(&[("a;size=3", "ACGT"), ("b;size=2", "ACGT"), ("c", "TTTT")]);
    let expected = fasta(&[("a;size=3", "ACGT"), ("c", "TTTT")]);
    assert_eq!(derep(&["--no-sizeout"], &input), expected);
}

#[test]
fn reads_phred64_and_writes_phred33() {
    let input = fastq(&[("a", "ACGT", "h@Jh"), ("b", "ACGT", "hhJ`")]);
    let expected = fastq(&[("a;size=2", "ACGT", "I%+C")]);
    assert_eq!(derep(&["--ascii", "64"], &input), expected);
    let expected = fastq(&[("a;size=2", "ACGT", "hDJb")]);
    assert_eq!(
        derep(&["--ascii", "64", "--ascii-out", "64"], &input),
        expected
    );
}

#[test]
fn limits_written_qualities() {
    let input = fastq(&[("a", "ACGTA", "!+5?I"), ("b", "ACGTA", "5++5~")]);
    let expected = fastq(&[("a;size=2", "ACGTA", "&+-7?")]);
    let args = ["--qmin-out", "5", "--qmax-out", "30"];
    assert_eq!(derep(&args, &input), expected);
}

#[test]
fn keeps_the_best_quality() {
    let input = fastq(&[("a", "ACGTA", "!+5?I"), ("b", "ACGTA", "5++5~")]);
    let expected = fastq(&[("a;size=2", "ACGTA", "5+5?~")]);
    assert_eq!(derep(&["--best-qual"], &input), expected);
}

#[test]
fn rejects_qualities_out_of_range() {
    let input = fastq(&[("a", "ACGT", "II5I")]);
    for (args, message) in [
        (
            &["--qmax", "30"][..],
            "quality score 40 above the maximum 30 in read a",
        ),
        (
            &["--qmin", "25"],
            "quality score 20 below the minimum 25 in read a",
        ),
        (
            &["--ascii", "64"],
            "quality score -11 below the minimum 0 in read a",
        ),
    ] {
        let output = run(args, &input);
        assert!(!output.status.success(), "{args:?} should fail");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(message),
            "{args:?}"
        );
    }
}

#[test]
fn rejects_invalid_quality_ranges() {
    let input = fastq(&[("a", "ACGT", "IIII")]);
    for args in [
        &["--qmin", "10", "--qmax", "5"][..],
        &["--ascii", "64", "--qmin", "-32"],
        &["--ascii", "64", "--qmax", "63"],
        &["--qmin-out", "10", "--qmax-out", "5"],
        &["--ascii-out", "64", "--qmax-out", "63"],
        &["--ascii", "40"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn rejects_quality_options_with_fasta() {
    let input = fasta(&[("a", "ACGT")]);
    for args in [
        &["--ascii", "64"][..],
        &["--qmin", "0"],
        &["--qmax-out", "41"],
        &["--best-qual"],
    ] {
        assert!(!run(args, &input).status.success(), "{args:?} should fail");
    }
}

#[test]
fn writes_uc_and_tables_of_every_cluster() {
    let (uc, table) = (temp_path("clusters.uc"), temp_path("clusters.tsv"));
    let input = fasta(&[
        ("a;size=2", "AACC"),
        ("b", "GGTT"),
        ("c", "TTTT"),
        ("d", "aacc"),
    ]);
    let output = derep(
        &[
            "--strand",
            "--uc",
            uc.to_str().unwrap(),
            "--tabbedout",
            table.to_str().unwrap(),
        ],
        &input,
    );
    let (uc_text, table_text) = (
        std::fs::read_to_string(&uc).unwrap(),
        std::fs::read_to_string(&table).unwrap(),
    );
    std::fs::remove_file(uc).unwrap();
    std::fs::remove_file(table).unwrap();
    assert_eq!(output, fasta(&[("a;size=4", "AACC"), ("c;size=1", "TTTT")]));
    assert_eq!(
        uc_text,
        "S\t0\t4\t*\t*\t*\t*\t*\ta;size=2\t*\n\
         H\t0\t4\t100.0\t-\t0\t0\t*\tb\ta;size=2\n\
         H\t0\t4\t100.0\t+\t0\t0\t*\td\ta;size=2\n\
         S\t1\t4\t*\t*\t*\t*\t*\tc\t*\n\
         C\t0\t4\t*\t*\t*\t*\t*\ta;size=2\t*\n\
         C\t1\t1\t*\t*\t*\t*\t*\tc\t*\n"
    );
    // the second column is the label of the cluster's first sequence
    assert_eq!(
        table_text,
        "a;size=2\ta;size=2\t0\t0\t3\ta;size=2\n\
         b\ta;size=2\t0\t1\t3\ta;size=2\n\
         d\ta;size=2\t0\t2\t3\ta;size=2\n\
         c\tc\t1\t0\t1\tc\n"
    );
}

#[test]
fn writes_whole_headers_in_uc_and_tables() {
    let (uc, table) = (temp_path("blanks.uc"), temp_path("blanks.tsv"));
    let input = fasta(&[("a d1", "AACC"), ("b d2", "AACC")]);
    let output = derep(
        &[
            "--uc",
            uc.to_str().unwrap(),
            "--tabbedout",
            table.to_str().unwrap(),
        ],
        &input,
    );
    let (uc_text, table_text) = (
        std::fs::read_to_string(&uc).unwrap(),
        std::fs::read_to_string(&table).unwrap(),
    );
    std::fs::remove_file(uc).unwrap();
    std::fs::remove_file(table).unwrap();
    assert_eq!(output, fasta(&[("a d1;size=2", "AACC")]));
    assert_eq!(
        uc_text,
        "S\t0\t4\t*\t*\t*\t*\t*\ta d1\t*\n\
         H\t0\t4\t100.0\t+\t0\t0\t*\tb d2\ta d1\n\
         C\t0\t2\t*\t*\t*\t*\t*\ta d1\t*\n"
    );
    assert_eq!(
        table_text,
        "a d1\ta d1\t0\t0\t2\ta d1\nb d2\ta d1\t0\t1\t2\ta d1\n"
    );
}

#[test]
fn merges_only_identical_labels_with_same_label() {
    // vsearch --derep_id: whole headers are compared
    let input = fasta(&[("a", "ACGT"), ("b", "ACGT"), ("a", "acgt"), ("a", "TTTT")]);
    let expected = fasta(&[
        ("a;size=2", "ACGT"),
        ("a;size=1", "TTTT"),
        ("b;size=1", "ACGT"),
    ]);
    assert_eq!(derep(&["--same-label"], &input), expected);
    let input = fasta(&[("a x", "ACGT"), ("a y", "ACGT"), ("a x", "acgt")]);
    let expected = fasta(&[("a x;size=2", "ACGT"), ("a y;size=1", "ACGT")]);
    assert_eq!(derep(&["--same-label"], &input), expected);
}

#[test]
fn merges_prefixes_into_longer_sequences() {
    // vsearch --derep_prefix
    let uc = temp_path("prefix.uc");
    let input = fasta(&[
        ("a", "ACGTAC"),
        ("b", "ACG"),
        ("c", "acgt"),
        ("d", "TTT"),
        ("e", "ACGTAC"),
    ]);
    let output = derep(&["--prefix", "--uc", uc.to_str().unwrap()], &input);
    let uc_text = std::fs::read_to_string(&uc).unwrap();
    std::fs::remove_file(uc).unwrap();
    assert_eq!(
        output,
        fasta(&[("a;size=4", "ACGTAC"), ("d;size=1", "TTT")])
    );
    // members keep their own lengths
    assert_eq!(
        uc_text,
        "S\t0\t6\t*\t*\t*\t*\t*\ta\t*\n\
         H\t0\t4\t100.0\t+\t0\t0\t*\tc\ta\n\
         H\t0\t3\t100.0\t+\t0\t0\t*\tb\ta\n\
         H\t0\t6\t100.0\t+\t0\t0\t*\te\ta\n\
         S\t1\t3\t*\t*\t*\t*\t*\td\t*\n\
         C\t0\t4\t*\t*\t*\t*\t*\ta\t*\n\
         C\t1\t1\t*\t*\t*\t*\t*\td\t*\n"
    );
    // an empty sequence is a prefix of every other
    let input = fasta(&[("a", "AAAA"), ("e", ""), ("b;size=3", "CCCC"), ("f", "")]);
    let expected = fasta(&[("b;size=5", "CCCC"), ("a;size=1", "AAAA")]);
    assert_eq!(derep(&["--prefix"], &input), expected);
}

#[test]
fn refuses_fastq_with_prefix() {
    // vsearch --derep_prefix writes only FASTA, and changing the format is
    // conv's job
    let input = fastq(&[("a", "ACGT", "IIII"), ("b", "AC", "II")]);
    assert!(!run(&["--prefix"], &input).status.success());
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--unknown"][..],
        &["--prefix", "--strand"],
        &["--prefix", "--same-label"],
    ] {
        assert!(!run(args, "").status.success(), "{args:?} should fail");
    }
}

#[test]
fn leaves_other_jobs_to_other_tools() {
    // filtering (sqfilt), topn, relabelling (relab), the output format
    // (conv) and the map that --tabbedout covers
    for args in [
        &["--min", "2"][..],
        &["--max", "2"],
        &["--min-len", "4"],
        &["--max-len", "4"],
        &["--top", "2"],
        &["--relabel", "U"],
        &["--relabel-self"],
        &["--relabel-sha1"],
        &["--relabel-md5"],
        &["--relabel-keep"],
        &["--xsize"],
        &["--xee"],
        &["--lengthout"],
        &["--xlength"],
        &["--sample", "S"],
        &["--label-suffix", "x"],
        &["--fasta"],
        &["--fasta-width", "80"],
        &["--map", "map.tsv"],
        &["--full-labels"],
    ] {
        let output = run(args, &fasta(&[("a", "ACGT")]));
        assert!(!output.status.success(), "{args:?} should be refused");
    }
}
#[test]
fn accepts_empty_input() {
    assert_eq!(derep(&[], ""), "");
}

/// Eighty reads cut from three templates at many lengths, so that most
/// are prefixes of others: some in lower case or with U for T, some with
/// an ambiguity code, some changed near their end (a branch), some empty,
/// with abundances.
fn prefix_reads() -> String {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = move |n: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % n
    };
    let templates: Vec<Vec<u8>> = (0..3)
        .map(|_| (0..40).map(|_| b"ACGT"[next(4) as usize]).collect())
        .collect();
    let mut records = Vec::new();
    for i in 0..80 {
        let template = &templates[next(3) as usize];
        let len = next(41) as usize;
        let mut seq = template[..len].to_vec();
        match next(8) {
            0 => seq.make_ascii_lowercase(),
            1 => seq
                .iter_mut()
                .filter(|c| **c == b'T')
                .for_each(|c| *c = b'U'),
            2 if len > 0 => seq[next(len as u64) as usize] = b'R',
            3 if len > 0 => seq[len - 1] = b"ACGT"[next(4) as usize],
            _ => {}
        }
        let label = match next(4) {
            0 => format!("r{i};size={}", 1 + next(5)),
            _ => format!("r{i}"),
        };
        records.push(format!(">{label}\n{}\n", String::from_utf8(seq).unwrap()));
    }
    records.concat()
}

#[test]
fn merges_many_prefixes_as_vsearch() {
    // vsearch --derep_prefix --sizein --sizeout --minseqlength 0
    // --notrunclabels --fasta_width 0 --output --uc
    let input = prefix_reads();
    let uc = temp_path("many-prefixes.uc");
    let output = derep(&["--prefix", "--uc", uc.to_str().unwrap()], &input);
    let uc_text = std::fs::read_to_string(&uc).unwrap();
    std::fs::remove_file(uc).unwrap();
    assert_eq!(output, MANY_PREFIXES);
    assert_eq!(uc_text, MANY_PREFIXES_UC);
}

const MANY_PREFIXES: &str = r##">r37;size=35
TATGGCACTCACCGATGACGGTTTAGTCTGCTGAGG
>r21;size=8
CACT
>r4;size=8
TATTCTGTGACTTGT
>r12;size=7
CACGGGAATATACAATGTTCR
>r38;size=7
CACGGGAATATACAATGTTCCCTATTACTAGCGCTAA
>r72;size=7
TATTCTGTGACTTGCATGCGCA
>r15;size=6
TATTCR
>r16;size=5
TATTCTGTGACTTGCATGCGCCAATAGTCGCAGAGG
>r5;size=3
CACGGGT
>r59;size=3
TATTCTGTGR
>r1;size=1
CACGGGAATATACAATGTTCCCTATTRCTAG
>r18;size=1
CACGGGAATATACAATGRTCCCTATTACT
>r23;size=1
TRTGGCACTCACC
>r34;size=1
CACGGGAATATARAATGTTCCCT
>r50;size=1
RATTCTGTGA
>r64;size=1
CACGGGAATATACARTGTTCCCTATTACTAGCGC
>r65;size=1
CACGGGAARATACAATGTTCCCTATTACTAG
"##;
const MANY_PREFIXES_UC: &str = r##"S	0	36	*	*	*	*	*	r37	*
H	0	35	100.0	+	0	0	*	r39	r37
H	0	34	100.0	+	0	0	*	r45	r37
H	0	33	100.0	+	0	0	*	r48	r37
H	0	30	100.0	+	0	0	*	r47	r37
H	0	27	100.0	+	0	0	*	r62	r37
H	0	23	100.0	+	0	0	*	r30	r37
H	0	22	100.0	+	0	0	*	r19;size=3	r37
H	0	19	100.0	+	0	0	*	r77;size=2	r37
H	0	18	100.0	+	0	0	*	r66	r37
H	0	16	100.0	+	0	0	*	r25;size=3	r37
H	0	15	100.0	+	0	0	*	r41	r37
H	0	14	100.0	+	0	0	*	r52	r37
H	0	13	100.0	+	0	0	*	r10	r37
H	0	12	100.0	+	0	0	*	r29;size=1	r37
H	0	11	100.0	+	0	0	*	r22;size=3	r37
H	0	8	100.0	+	0	0	*	r11	r37
H	0	5	100.0	+	0	0	*	r69;size=3	r37
H	0	13	100.0	+	0	0	*	r6	r37
H	0	19	100.0	+	0	0	*	r17;size=1	r37
H	0	22	100.0	+	0	0	*	r61	r37
H	0	23	100.0	+	0	0	*	r58	r37
H	0	23	100.0	+	0	0	*	r78	r37
H	0	35	100.0	+	0	0	*	r42	r37
H	0	35	100.0	+	0	0	*	r63	r37
H	0	36	100.0	+	0	0	*	r73	r37
S	1	4	*	*	*	*	*	r21;size=3	*
H	1	3	100.0	+	0	0	*	r51	r21;size=3
H	1	1	100.0	+	0	0	*	r57	r21;size=3
H	1	0	100.0	+	0	0	*	r24;size=1	r21;size=3
H	1	0	100.0	+	0	0	*	r8	r21;size=3
H	1	3	100.0	+	0	0	*	r60	r21;size=3
S	2	15	*	*	*	*	*	r4	*
H	2	14	100.0	+	0	0	*	r31;size=3	r4
H	2	13	100.0	+	0	0	*	r71;size=1	r4
H	2	11	100.0	+	0	0	*	r20	r4
H	2	10	100.0	+	0	0	*	r76	r4
H	2	14	100.0	+	0	0	*	r56	r4
S	3	21	*	*	*	*	*	r12	*
H	3	20	100.0	+	0	0	*	r40	r12
H	3	19	100.0	+	0	0	*	r43	r12
H	3	14	100.0	+	0	0	*	r36	r12
H	3	13	100.0	+	0	0	*	r0	r12
H	3	9	100.0	+	0	0	*	r54	r12
H	3	19	100.0	+	0	0	*	r75	r12
S	4	37	*	*	*	*	*	r38	*
H	4	32	100.0	+	0	0	*	r55	r38
H	4	31	100.0	+	0	0	*	r26;size=1	r38
H	4	30	100.0	+	0	0	*	r27	r38
H	4	29	100.0	+	0	0	*	r53	r38
H	4	27	100.0	+	0	0	*	r35	r38
H	4	22	100.0	+	0	0	*	r70	r38
S	5	22	*	*	*	*	*	r72	*
H	5	19	100.0	+	0	0	*	r44	r72
H	5	18	100.0	+	0	0	*	r7;size=4	r72
H	5	16	100.0	+	0	0	*	r33	r72
S	6	6	*	*	*	*	*	r15	*
H	6	5	100.0	+	0	0	*	r79	r15
H	6	4	100.0	+	0	0	*	r9	r15
H	6	3	100.0	+	0	0	*	r2	r15
H	6	2	100.0	+	0	0	*	r28	r15
H	6	1	100.0	+	0	0	*	r68	r15
S	7	36	*	*	*	*	*	r16	*
H	7	30	100.0	+	0	0	*	r14	r16
H	7	29	100.0	+	0	0	*	r67	r16
H	7	26	100.0	+	0	0	*	r32	r16
H	7	23	100.0	+	0	0	*	r46	r16
S	8	7	*	*	*	*	*	r5	*
H	8	4	100.0	+	0	0	*	r13	r5
H	8	4	100.0	+	0	0	*	r49	r5
S	9	10	*	*	*	*	*	r59	*
H	9	7	100.0	+	0	0	*	r3	r59
H	9	6	100.0	+	0	0	*	r74	r59
S	10	31	*	*	*	*	*	r1	*
S	11	29	*	*	*	*	*	r18	*
S	12	13	*	*	*	*	*	r23	*
S	13	23	*	*	*	*	*	r34	*
S	14	10	*	*	*	*	*	r50	*
S	15	34	*	*	*	*	*	r64	*
S	16	31	*	*	*	*	*	r65	*
C	0	35	*	*	*	*	*	r37	*
C	1	8	*	*	*	*	*	r21;size=3	*
C	2	8	*	*	*	*	*	r4	*
C	3	7	*	*	*	*	*	r12	*
C	4	7	*	*	*	*	*	r38	*
C	5	7	*	*	*	*	*	r72	*
C	6	6	*	*	*	*	*	r15	*
C	7	5	*	*	*	*	*	r16	*
C	8	3	*	*	*	*	*	r5	*
C	9	3	*	*	*	*	*	r59	*
C	10	1	*	*	*	*	*	r1	*
C	11	1	*	*	*	*	*	r18	*
C	12	1	*	*	*	*	*	r23	*
C	13	1	*	*	*	*	*	r34	*
C	14	1	*	*	*	*	*	r50	*
C	15	1	*	*	*	*	*	r64	*
C	16	1	*	*	*	*	*	r65	*
"##;
