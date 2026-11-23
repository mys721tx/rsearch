//! End-to-end tests of the fqstat binary.
//!
//! Expected reports are those `vsearch --fastq_stats` writes to its log
//! (`--log`), without the lines the log adds before and after. They give
//! the distribution of read lengths and of quality scores, and the average
//! quality and expected errors at each length. They also give the reads
//! that would pass four expected-error cutoffs and four quality cutoffs at
//! each length, and the totals. `--ascii`, `--qmin` and `--qmax` are
//! vsearch's `--fastq_ascii`, `--fastq_qmin` and `--fastq_qmax`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// Runs fqstat with `args`, feeding `input` on stdin.
fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_fqstat"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("fqstat should start");
    let mut stdin = child.stdin.take().unwrap();
    let fed = input.as_bytes().to_vec();
    let feeder = std::thread::spawn(move || {
        let _ = stdin.write_all(&fed);
    });
    let output = child.wait_with_output().unwrap();
    feeder.join().unwrap();
    output
}

/// Runs fqstat and returns its report, which must succeed.
fn fqstat(args: &[&str], input: &str) -> String {
    let output = run(args, input);
    assert!(
        output.status.success(),
        "fqstat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
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
    std::env::temp_dir().join(format!("rsearch-fqstat-{}-{name}", std::process::id()))
}

/// The headers of the report's sections, each followed by its rows.
const LENGTHS: &str = "
Read length distribution
      L           N      Pct   AccPct
-------  ----------  -------  -------
";
const SCORES: &str = "
Q score distribution
ASCII    Q       Pe           N      Pct   AccPct
-----  ---  -------  ----------  -------  -------
";
const AVERAGES: &str = "
    L  PctRecs  AvgQ  P(AvgQ)      AvgP  AvgEE       Rate   RatePct
-----  -------  ----  -------  --------  -----  ---------  --------
";
const EE_CUTOFFS: &str = "
    L   1.0000   0.5000   0.2500   0.1000   1.0000   0.5000   0.2500   0.1000
-----  -------  -------  -------  -------  -------  -------  -------  -------
";
const Q_CUTOFFS: &str = "
Truncate at first Q
  Len     Q=5    Q=10    Q=15    Q=20
-----  ------  ------  ------  ------
";

/// Assembles a report from the rows of its sections and its totals.
fn report(sections: [&str; 5], totals: &str) -> String {
    let headers = [LENGTHS, SCORES, AVERAGES, EE_CUTOFFS, Q_CUTOFFS];
    let body: String = headers
        .iter()
        .zip(sections)
        .map(|(header, rows)| format!("{header}{rows}"))
        .collect();
    format!("{body}\n{totals}")
}

#[test]
fn reports_lengths_scores_and_cutoffs() {
    let input = fastq(&[
        ("a", "ACGTACGTAC", "IIIII5555#"),
        ("b", "ACGTAC", "III+++"),
        ("c", "AC", "!I"),
        ("d", "ACGTACGTAC", "IIIIIIIIII"),
    ]);
    let expected = report(
        [
            "\
>=   10           2    50.0%    50.0%
      6           1    25.0%    75.0%
      2           1    25.0%   100.0%
",
            "    I   40  0.00010          19    67.9%    67.9%
    5   20  0.01000           4    14.3%    82.1%
    +   10  0.10000           3    10.7%    92.9%
    #    2  0.63096           1     3.6%    96.4%
    !    0  1.00000           1     3.6%   100.0%
",
            "    2   100.0%  40.0  0.00010  0.000100   0.25   0.125087   12.509%
    3    75.0%  40.0  0.00010  0.000100   0.00   0.000100    0.010%
    4    75.0%  30.0  0.00100  0.033400   0.03   0.008425    0.843%
    5    75.0%  30.0  0.00100  0.033400   0.07   0.013420    1.342%
    6    75.0%  23.3  0.00464  0.036700   0.10   0.017300    1.730%
    7    50.0%  30.0  0.00100  0.005050   0.01   0.001514    0.151%
    8    50.0%  30.0  0.00100  0.005050   0.02   0.001956    0.196%
    9    50.0%  30.0  0.00100  0.005050   0.02   0.002300    0.230%
   10    50.0%  21.0  0.00794  0.315529   0.34   0.033623    3.362%
",
            "   10        2        1        1        1   50.00%   25.00%   25.00%   25.00%
    9        2        2        2        2   50.00%   50.00%   50.00%   50.00%
    8        2        2        2        2   50.00%   50.00%   50.00%   50.00%
    7        2        2        2        2   50.00%   50.00%   50.00%   50.00%
    6        3        3        2        2   75.00%   75.00%   50.00%   50.00%
    5        3        3        3        2   75.00%   75.00%   75.00%   50.00%
    4        3        3        3        2   75.00%   75.00%   75.00%   50.00%
    3        3        3        3        3   75.00%   75.00%   75.00%   75.00%
    2        3        3        3        3   75.00%   75.00%   75.00%   75.00%
    1        4        3        3        3  100.00%   75.00%   75.00%   75.00%
",
            "   10   25.0%   25.0%   25.0%   25.0%
    9   50.0%   50.0%   50.0%   25.0%
    8   50.0%   50.0%   50.0%   25.0%
    7   50.0%   50.0%   50.0%   25.0%
    6   75.0%   50.0%   50.0%   25.0%
    5   75.0%   50.0%   50.0%   50.0%
",
        ],
        "         4  Recs (0.0M), 0 too long
       7.0  Avg length
      0.0M  Bases
",
    );
    assert_eq!(fqstat(&[], &input), expected);
}

#[test]
fn stops_the_expected_error_table_where_no_read_passes() {
    // the expected-error table starts at the last length where a read still
    // has at most 1.0 expected errors. Ten scores of 10 sum to just under
    // 1.0. The quality table covers the longer half of the lengths
    let input = fastq(&[
        ("a", "ACGTACGTACGTAC", "++++++++++++++"),
        ("b", "ACGTACGTACGT", "5555555555++"),
        ("c", "ACGTACGTACGTACGTACGTA", "+++++IIIII+++++IIIII#"),
    ]);
    let expected = report(
        [
            "\
>=   21           1    33.3%    33.3%
     14           1    33.3%    66.7%
     12           1    33.3%   100.0%
",
            "    I   40  0.00010          10    21.3%    21.3%
    5   20  0.01000          10    21.3%    42.6%
    +   10  0.10000          26    55.3%    97.9%
    #    2  0.63096           1     2.1%   100.0%
",
            "    2   100.0%  13.3  0.04642  0.070000   0.14   0.070000    7.000%
    3   100.0%  13.3  0.04642  0.070000   0.21   0.070000    7.000%
    4   100.0%  13.3  0.04642  0.070000   0.28   0.070000    7.000%
    5   100.0%  13.3  0.04642  0.070000   0.35   0.070000    7.000%
    6   100.0%  23.3  0.00464  0.036700   0.39   0.064450    6.445%
    7   100.0%  23.3  0.00464  0.036700   0.42   0.060486    6.049%
    8   100.0%  23.3  0.00464  0.036700   0.46   0.057512    5.751%
    9   100.0%  23.3  0.00464  0.036700   0.50   0.055200    5.520%
   10   100.0%  23.3  0.00464  0.036700   0.53   0.053350    5.335%
   11   100.0%  10.0  0.10000  0.100000   0.63   0.057591    5.759%
   12   100.0%  10.0  0.10000  0.100000   0.73   0.061125    6.112%
   13    66.7%  10.0  0.10000  0.100000   1.05   0.080788    8.079%
   14    66.7%  10.0  0.10000  0.100000   1.15   0.082161    8.216%
   15    33.3%  10.0  0.10000  0.100000   1.00   0.066700    6.670%
   16    33.3%  40.0  0.00010  0.000100   1.00   0.062537    6.254%
   17    33.3%  40.0  0.00010  0.000100   1.00   0.058865    5.886%
   18    33.3%  40.0  0.00010  0.000100   1.00   0.055600    5.560%
   19    33.3%  40.0  0.00010  0.000100   1.00   0.052679    5.268%
   20    33.3%  40.0  0.00010  0.000100   1.00   0.050050    5.005%
   21    33.3%   2.0  0.63096  0.630957   1.63   0.077712    7.771%
",
            "   14        1        0        0        0   33.33%    0.00%    0.00%    0.00%
   13        1        0        0        0   33.33%    0.00%    0.00%    0.00%
   12        2        1        0        0   66.67%   33.33%    0.00%    0.00%
   11        2        1        1        0   66.67%   33.33%   33.33%    0.00%
   10        3        1        1        1  100.00%   33.33%   33.33%   33.33%
    9        3        1        1        1  100.00%   33.33%   33.33%   33.33%
    8        3        1        1        1  100.00%   33.33%   33.33%   33.33%
    7        3        1        1        1  100.00%   33.33%   33.33%   33.33%
    6        3        1        1        1  100.00%   33.33%   33.33%   33.33%
    5        3        3        1        1  100.00%  100.00%   33.33%   33.33%
    4        3        3        1        1  100.00%  100.00%   33.33%   33.33%
    3        3        3        1        1  100.00%  100.00%   33.33%   33.33%
    2        3        3        3        1  100.00%  100.00%  100.00%   33.33%
    1        3        3        3        3  100.00%  100.00%  100.00%  100.00%
",
            "   21    0.0%    0.0%    0.0%    0.0%
   20   33.3%    0.0%    0.0%    0.0%
   19   33.3%    0.0%    0.0%    0.0%
   18   33.3%    0.0%    0.0%    0.0%
   17   33.3%    0.0%    0.0%    0.0%
   16   33.3%    0.0%    0.0%    0.0%
   15   33.3%    0.0%    0.0%    0.0%
   14   66.7%    0.0%    0.0%    0.0%
   13   66.7%    0.0%    0.0%    0.0%
   12  100.0%    0.0%    0.0%    0.0%
   11  100.0%    0.0%    0.0%    0.0%
   10  100.0%   33.3%   33.3%    0.0%
",
        ],
        "         3  Recs (0.0M), 0 too long
      15.7  Avg length
      0.0M  Bases
",
    );
    assert_eq!(fqstat(&[], &input), expected);
}

#[test]
fn reports_empty_input() {
    let expected = report(
        ["", "", "", "", ""],
        "         0  Recs (0.0M), 0 too long
      0.0M  Bases
",
    );
    assert_eq!(fqstat(&[], ""), expected);
}

#[test]
fn reports_empty_reads() {
    let expected = report(
        [">=    0           1   100.0%   100.0%\n", "", "", "", ""],
        "         1  Recs (0.0M), 0 too long
       0.0  Avg length
      0.0M  Bases
",
    );
    assert_eq!(fqstat(&[], &fastq(&[("e", "", "")])), expected);
    let expected = report(
        [
            ">=    1           1    50.0%    50.0%
      0           1    50.0%   100.0%
",
            "    I   40  0.00010           1   100.0%   100.0%\n",
            "",
            "    1        1        1        1        1   50.00%   50.00%   50.00%   50.00%\n",
            "    1   50.0%   50.0%   50.0%   50.0%\n",
        ],
        "         2  Recs (0.0M), 0 too long
       0.5  Avg length
      0.0M  Bases
",
    );
    let input = fastq(&[("e", "", ""), ("a", "A", "I")]);
    assert_eq!(fqstat(&[], &input), expected);
}

#[test]
fn reads_other_offsets_and_negative_scores() {
    let expected = report(
        [
            ">=    4           1   100.0%   100.0%\n",
            "    h   40  0.00010           2    50.0%    50.0%
    @    0  1.00000           1    25.0%    75.0%
    ;   -5  3.16228           1    25.0%   100.0%
",
            "    2   100.0%   0.0  1.00000  1.000000   4.16   2.081139  208.114%
    3   100.0%  40.0  0.00010  0.000100   4.16   1.387459  138.746%
    4   100.0%  40.0  0.00010  0.000100   4.16   1.040619  104.062%
",
            "",
            "    4    0.0%    0.0%    0.0%    0.0%
    3    0.0%    0.0%    0.0%    0.0%
    2    0.0%    0.0%    0.0%    0.0%
",
        ],
        "         1  Recs (0.0M), 0 too long
       4.0  Avg length
      0.0M  Bases
",
    );
    let input = fastq(&[("a", "ACGT", ";@hh")]);
    assert_eq!(fqstat(&["--ascii", "64", "--qmin", "-5"], &input), expected);
}

#[test]
fn checks_the_range_of_scores() {
    // scores up to 126 - offset are accepted by default
    let input = fastq(&[("a", "AC", "IK")]);
    assert!(fqstat(&[], &input).contains("    K   42  0.00006"));
    // '!' is -31 at offset 64
    let low = fastq(&[("a", "AC", "I!")]);
    for (args, input, message) in [
        (&["--qmax", "41"][..], &input, "above"),
        (&["--qmin", "41"], &input, "below"),
        (&["--ascii", "64"], &low, "below"),
    ] {
        let output = run(args, input);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains(message), "{args:?}: {stderr}");
    }
}

#[test]
fn ignores_abundances() {
    // vsearch --fastq_stats reads no abundance and checks none
    let output = run(&[], &fastq(&[("a;size=0", "AC", "II")]));
    assert!(output.status.success());
}

#[test]
fn rejects_fasta() {
    let output = run(&[], ">a\nACGT\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("FASTQ file expected"), "{stderr}");
}

#[test]
fn rejects_invalid_input() {
    let output = run(&[], "@a\nACGT\nIIII\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let output = run(&[], "@a\nAC\n+\nI\x1f\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn rejects_invalid_arguments() {
    for args in [
        &["--unknown"][..],
        &["--ascii", "50"],
        &["--qmin", "10", "--qmax", "5"],
        &["--qmin", "-1"],
        &["--ascii", "64", "--qmax", "63"],
    ] {
        let output = run(args, &fastq(&[("a", "AC", "II")]));
        assert_eq!(output.status.code(), Some(2), "{args:?} should fail");
    }
}

#[test]
fn reads_and_writes_files() {
    let (input, output) = (temp_path("in.fq"), temp_path("out.txt"));
    std::fs::write(&input, fastq(&[("a", "AC", "II")])).unwrap();
    let stdout = fqstat(
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
        ],
        "",
    );
    let written = std::fs::read_to_string(&output).unwrap();
    std::fs::remove_file(input).unwrap();
    std::fs::remove_file(output).unwrap();
    assert_eq!(stdout, "");
    assert!(written.starts_with(LENGTHS), "{written}");
    assert!(written.ends_with(
        "         1  Recs (0.0M), 0 too long\n       2.0  Avg length\n      0.0M  Bases\n"
    ));
}
