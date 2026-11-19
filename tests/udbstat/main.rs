//! End-to-end tests of the udbstat binary.
//!
//! Expected reports are what vsearch `--udbstats` writes to `--log` for the
//! same UDB files, written here by mkudb. A report gives the summary of the
//! database and the word statistics. It lists the eleven words held by the
//! most sequences with the first eight of them, and the distribution of the
//! words over sequence counts. udbstat writes the report to `-o`. It reads
//! the whole file with vsearch's checks, and reads stdin too, where vsearch
//! needs a file.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-udbstat-{}-{name}", std::process::id()))
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

/// The UDB file mkudb writes for `fasta` at word length `word_len`.
fn mkudb(word_len: u32, fasta: &str) -> Vec<u8> {
    let word_len = word_len.to_string();
    let output = run(
        env!("CARGO_BIN_EXE_mkudb"),
        &["--word-len", &word_len],
        fasta.as_bytes(),
    );
    assert!(output.status.success());
    output.stdout
}

/// Runs udbstat on `input`.
fn udbstat(input: &[u8]) -> Output {
    run(env!("CARGO_BIN_EXE_udbstat"), &[], input)
}

/// Runs udbstat, which must succeed, and returns its report.
fn report(input: &[u8]) -> String {
    let output = udbstat(input);
    assert!(
        output.status.success(),
        "udbstat failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    String::from_utf8(output.stdout).unwrap()
}

/// Twenty sequences of ten bases sharing their first six.
const SHARED: [&str; 20] = [
    "ACGTACTTTC",
    "ACGTACCTCA",
    "ACGTACTGCA",
    "ACGTACATTC",
    "ACGTACAAAA",
    "ACGTACCCAT",
    "ACGTACGTCC",
    "ACGTACGTAA",
    "ACGTACTGTA",
    "ACGTACGGCG",
    "ACGTACAAAT",
    "ACGTACAGTA",
    "ACGTACAACC",
    "ACGTACATTT",
    "ACGTACTACG",
    "ACGTACGAGG",
    "ACGTACATAC",
    "ACGTACCAAA",
    "ACGTACTTCC",
    "ACGTACTCCT",
];

fn shared() -> String {
    SHARED
        .iter()
        .enumerate()
        .map(|(i, s)| format!(">r{i}\n{s}\n"))
        .collect()
}

#[test]
fn reports_shared_words() {
    // eight sequence numbers at most, then an ellipsis; the median of the
    // counts; the buckets up to the number of sequences
    assert_eq!(
        report(&mkudb(3, &shared())),
        r#"200 nt in 20 seqs, min 10, max 10, avg 10

      Alphabet  nt
    Word width  3
     Word ones  3
        Spaced  No
        Hashed  No
         Coded  No
       Stepped  No
         Slots  64 (0.1k)
       DBAccel  100%

       200  DB size (0.2k)
       147  Words
         1  Median size
       2.3  Mean size

     iWord         sWord         Cap        Size  Row
----------  ------------  ----------  ----------  ---
         6           ACG           0          20  0 1 2 3 4 5 6 7...
        27           CGT           0          20  0 1 2 3 4 5 6 7...
        44           GTA           0          20  0 1 2 3 4 5 6 7...
        49           TAC           0          20  0 1 2 3 4 5 6 7...
         4           ACA           0           7  3 4 10 11 12 13 16
         7           ACT           0           6  0 2 8 14 18 19
         5           ACC           0           4  1 5 12 17
        16           CAA           0           4  4 10 12 17
        19           CAT           0           4  3 5 13 16
         0           AAA           0           3  4 10 17
        53           TCC           0           3  6 18 19


Word width  3
Slots       64
Words       147
Max size    20 (ACG)

   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct
----------  ----------  ----------  ----------  ------  ------
                     0         0.0        26.0   40.6%   40.6%  **************
                     1        19.0        19.0   29.7%   70.3%  **********
                     2        14.0         7.0   10.9%   81.2%  ****
         3           4        21.0         6.0    9.4%   90.6%  ***
         5           8        13.0         2.0    3.1%   93.8%  *
         9          16         0.0         0.0    0.0%   93.8%
        17          20        80.0         4.0    6.2%  100.0%  **
----------  ----------  ----------  ----------
                             147.0        64.0

       200  Upper
         0  Lower (0.0%)
       200  Total
       147  Indexed words
"#
    );
}

#[test]
fn reports_the_default_word_length() {
    // words never found ranked by their code; 65,536 words in thousands
    assert_eq!(
        report(&mkudb(8, &shared())),
        r#"200 nt in 20 seqs, min 10, max 10, avg 10

      Alphabet  nt
    Word width  8
     Word ones  8
        Spaced  No
        Hashed  No
         Coded  No
       Stepped  No
         Slots  65536 (65.5k)
       DBAccel  100%

       200  DB size (0.2k)
        60  Words
         0  Median size
       0.0  Mean size

     iWord         sWord         Cap        Size  Row
----------  ------------  ----------  ----------  ---
      6928      ACGTACAA           0           3  4 10 12
      6931      ACGTACAT           0           3  3 13 16
      6939      ACGTACGT           0           2  6 7
      6942      ACGTACTG           0           2  2 8
      6943      ACGTACTT           0           2  0 18
     27712      CGTACAAA           0           2  4 10
     27727      CGTACATT           0           2  3 13
      6930      ACGTACAG           0           1  11
      6932      ACGTACCA           0           1  17
      6933      ACGTACCC           0           1  5
      6935      ACGTACCT           0           1  1


Word width  8
Slots       65536
Words       60
Max size    3 (ACGTACAA)

   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct
----------  ----------  ----------  ----------  ------  ------
                     0         0.0       65.5k   99.9%   99.9%  *********************************
                     1        44.0        44.0    0.1%  100.0%
                     2        10.0         5.0    0.0%  100.0%
         3           4         6.0         2.0    0.0%  100.0%
         5           8         0.0         0.0    0.0%  100.0%
         9          16         0.0         0.0    0.0%  100.0%
        17          20         0.0         0.0    0.0%  100.0%
----------  ----------  ----------  ----------
                              60.0       65.5k

       200  Upper
         0  Lower (0.0%)
       200  Total
        60  Indexed words
"#
    );
}

#[test]
fn reports_one_sequence() {
    // the buckets stop below the number of sequences, and there is one
    // bucket only
    let input = mkudb(3, ">only\nACGTACGTAAC\n");
    assert_eq!(
        report(&input),
        r#"11 nt in 1 seqs, min 11, max 11, avg 11

      Alphabet  nt
    Word width  3
     Word ones  3
        Spaced  No
        Hashed  No
         Coded  No
       Stepped  No
         Slots  64 (0.1k)
       DBAccel  100%

        11  DB size (0.0k)
         6  Words
         0  Median size
       0.1  Mean size

     iWord         sWord         Cap        Size  Row
----------  ------------  ----------  ----------  ---
         1           AAC           0           1  0
         6           ACG           0           1  0
        27           CGT           0           1  0
        44           GTA           0           1  0
        48           TAA           0           1  0
        49           TAC           0           1  0
         0           AAA           0           0 
         2           AAG           0           0 
         3           AAT           0           0 
         4           ACA           0           0 
         5           ACC           0           0 


Word width  3
Slots       64
Words       6
Max size    1 (AAC)

   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct
----------  ----------  ----------  ----------  ------  ------
                     0         0.0        58.0   90.6%   90.6%  ******************************
----------  ----------  ----------  ----------
                               6.0        64.0

        11  Upper
         0  Lower (0.0%)
        11  Total
         6  Indexed words
"#
    );
}

/// Three hundred random sequences of 150 bases, from xorshift64.
fn random_sequences() -> String {
    let mut out = String::new();
    for i in 0..300u64 {
        let mut x = i.wrapping_mul(2_654_435_761).wrapping_add(1);
        let mut seq = String::new();
        for _ in 0..150 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            seq.push(b"ACGT"[(x >> 62) as usize] as char);
        }
        out.push_str(&format!(">q{i}\n{seq}\n"));
    }
    out
}

#[test]
fn reports_large_sizes_in_thousands() {
    assert_eq!(
        report(&mkudb(8, &random_sequences())),
        r#"45000 nt in 300 seqs, min 150, max 150, avg 150

      Alphabet  nt
    Word width  8
     Word ones  8
        Spaced  No
        Hashed  No
         Coded  No
       Stepped  No
         Slots  65536 (65.5k)
       DBAccel  100%

     45000  DB size (45.0k)
     42853  Words
         0  Median size
       0.7  Mean size

     iWord         sWord         Cap        Size  Row
----------  ------------  ----------  ----------  ---
     42431      GGCCGTTT           0           8  32 40 110 172 179 204 244 257
     11206      AGGTTACG           0           6  12 33 102 197 212 271
     12907      ATAGCGGT           0           6  14 19 64 209 217 298
     24072      CCTGAAGA           0           6  33 34 206 239 255 266
     48901      GTTTAACC           0           6  19 77 84 144 190 220
       178      AAAAGTAG           0           5  65 148 203 236 298
      2801      AAGGTTAC           0           5  33 120 197 271 287
      5737      ACCGCGGC           0           5  7 60 134 224 232
      5957      ACCTCACC           0           5  7 61 107 174 238
      6399      ACGATTTT           0           5  64 163 172 205 264
      7098      ACGTGTGG           0           5  126 136 149 285 292


Word width  8
Slots       65536
Words       42853
Max size    8 (GGCCGTTT)

   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct
----------  ----------  ----------  ----------  ------  ------
                     0         0.0       34.0k   52.0%   52.0%  *****************
                     1       22.4k       22.4k   34.1%   86.1%  ***********
                     2       14.5k      7238.0   11.0%   97.1%  ****
         3           4      5815.0      1855.0    2.8%   99.9%  *
         5           8       207.0        40.0    0.1%  100.0%
         9          16         0.0         0.0    0.0%  100.0%
        17          32         0.0         0.0    0.0%  100.0%
        33          64         0.0         0.0    0.0%  100.0%
        65         128         0.0         0.0    0.0%  100.0%
       129         256         0.0         0.0    0.0%  100.0%
       257         300         0.0         0.0    0.0%  100.0%
----------  ----------  ----------  ----------
                             42.9k       65.5k

     45000  Upper
         0  Lower (0.0%)
     45000  Total
     42853  Indexed words
"#
    );
}

/// The UDB file of one sequence at word length 3: the header (0), the
/// counts (200), UDB3 (456), six sequence numbers (460), the second header
/// (484), one header offset (516), "only\0" (520), one length (525) and
/// eleven residues (529), 540 bytes in all.
fn one() -> Vec<u8> {
    let udb = mkudb(3, ">only\nACGTACGTAAC\n");
    assert_eq!(udb.len(), 540);
    udb
}

/// `udb` with the 32-bit word at `offset` set to `value`.
fn with(udb: &[u8], offset: usize, value: u32) -> Vec<u8> {
    let mut udb = udb.to_vec();
    udb[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    udb
}

/// Runs udbstat on `input`, which must fail with `message`.
fn fails(input: &[u8], message: &str) {
    let output = udbstat(input);
    assert_eq!(output.status.code(), Some(1), "{message}");
    assert!(output.stdout.is_empty());
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        format!("Error: {message}\n")
    );
}

const SHORT: &str = "invalid UDB file: shorter than its header";
const TRUNCATED: &str = "invalid UDB file: truncated";
const INVALID: &str = "invalid UDB file";
const TRAILING: &str = "invalid UDB file: data after its end";

#[test]
fn rejects_truncated_files() {
    let udb = one();
    fails(b"", SHORT);
    fails(&udb[..199], SHORT);
    // in the counts, the sequence numbers and the sequences
    for end in [300, 470, 535, 539] {
        fails(&udb[..end], TRUNCATED);
    }
}

#[test]
fn rejects_data_after_the_end() {
    let mut udb = one();
    udb.push(b'x');
    fails(&udb, TRAILING);
}

#[test]
fn rejects_invalid_headers() {
    let udb = one();
    // no sequences, or more than the file holds 4 bytes for
    fails(&with(&udb, 52, 0), INVALID);
    fails(&with(&udb, 52, 135), INVALID);
    fails(&with(&udb, 52, 136), INVALID);
    // the second header: its magic numbers and the number of sequences
    for offset in [484, 488, 512] {
        fails(&with(&udb, offset, 0), INVALID);
    }
    fails(&with(&udb, 492, 2), INVALID);
}

#[test]
fn rejects_invalid_word_lists() {
    let udb = one();
    // a word held by more sequences than there are
    fails(&with(&udb, 200, 2), INVALID);
    // UDB3 missing
    fails(&with(&udb, 456, 0), INVALID);
    // a sequence number out of range
    fails(&with(&udb, 460, 1), INVALID);
    // more sequence numbers than the file could hold: two sequences at
    // word length 5, every word held by both
    let mut udb = mkudb(5, ">a\nACGTACGTAAC\n>b\nGGGTTTAC\n");
    for word in 0..1024 {
        udb = with(&udb, 200 + 4 * word, 2);
    }
    fails(&udb, INVALID);
}

#[test]
fn rejects_invalid_sequences() {
    let udb = one();
    // a header offset past the headers
    fails(&with(&udb, 516, 5), INVALID);
    // header bytes that do not end at the lengths
    fails(&with(&udb, 504, 4), INVALID);
    fails(&with(&udb, 504, 6), INVALID);
    // lengths that add up to more, or less, than the residues
    fails(&with(&udb, 525, 12), INVALID);
    fails(&with(&udb, 525, 10), INVALID);
    fails(&with(&udb, 496, 12), INVALID);
}

#[test]
fn reads_and_writes_files() {
    let (input, out) = (temp_path("in.udb"), temp_path("out.txt"));
    std::fs::write(&input, one()).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_udbstat"))
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&out);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(written, report(&one()));
}

#[test]
fn rejects_missing_files() {
    let output = Command::new(env!("CARGO_BIN_EXE_udbstat"))
        .args(["-i", temp_path("missing.udb").to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}
