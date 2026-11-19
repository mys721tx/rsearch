//! End-to-end tests of the udbinfo binary.
//!
//! Expected reports are those vsearch `--udbinfo` writes to `--log` (and to
//! stderr) for the same files. They hold the fields of the header of 50
//! words, as the file holds them. Like vsearch, udbinfo checks only the
//! header, and that the file is long enough for its number of sequences (4
//! bytes each). Unlike vsearch, which needs a file it can measure, it also
//! reads a UDB file from stdin.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-udbinfo-{}-{name}", std::process::id()))
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

/// Runs udbinfo on `input`.
fn udbinfo(input: &[u8]) -> Output {
    run(env!("CARGO_BIN_EXE_udbinfo"), &[], input)
}

/// Runs udbinfo, which must succeed, and returns its report.
fn report(input: &[u8]) -> String {
    let output = udbinfo(input);
    assert!(
        output.status.success(),
        "udbinfo failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    String::from_utf8(output.stdout).unwrap()
}

/// The header of a UDB file as mkudb writes it, for one sequence at word
/// length 8, with `fields` (word number, value) changed.
fn header(fields: &[(usize, u32)]) -> Vec<u8> {
    let mut words = [0u32; 50];
    words[0] = 0x5544_4246;
    words[2] = 32;
    words[4] = 8;
    words[5] = 1;
    words[6] = 100;
    words[13] = 1;
    words[17] = 0x0000_746e;
    words[49] = 0x5544_4266;
    for &(i, value) in fields {
        words[i] = value;
    }
    words.iter().flat_map(|w| w.to_le_bytes()).collect()
}

/// vsearch's report.
fn expected(seqs: u32, word_len: u32, slots: u32, step: u32, accel: u32, dict: &str) -> String {
    format!(
        "           Seqs  {seqs}\n     SeqIx bits  32\n          Alpha  nt (4)\n     Word width  \
         {word_len}\n          Slots  {slots}\n      Dict size  {dict}\n         DBstep  \
         {step}\n        DBAccel  {accel}%\n"
    )
}

#[test]
fn reports_mkudb_files() {
    let input = ">s1 desc\nACGTACGTACGT\n>s2\nAC\n>s3\nGGGGCCCCAAAATTTT\n";
    let udb = run(env!("CARGO_BIN_EXE_mkudb"), &[], input.as_bytes());
    assert!(udb.status.success());
    assert_eq!(
        report(&udb.stdout),
        expected(3, 8, 0, 1, 100, "65536 (65.5k)")
    );
    let udb = run(
        env!("CARGO_BIN_EXE_mkudb"),
        &["--word-len", "3"],
        input.as_bytes(),
    );
    assert_eq!(report(&udb.stdout), expected(3, 3, 0, 1, 100, "64 (0.1k)"));
}

#[test]
fn reports_dictionary_sizes() {
    // 4^k words, and in thousands with one decimal
    let sizes = [
        (3, "64 (0.1k)"),
        (4, "256 (0.3k)"),
        (5, "1024 (1.0k)"),
        (6, "4096 (4.1k)"),
        (7, "16384 (16.4k)"),
        (8, "65536 (65.5k)"),
        (9, "262144 (262.1k)"),
        (10, "1048576 (1048.6k)"),
        (11, "4194304 (4194.3k)"),
        (12, "16777216 (16777.2k)"),
        (13, "67108864 (67108.9k)"),
        (14, "268435456 (268435.5k)"),
        (15, "1073741824 (1073741.8k)"),
    ];
    for (word_len, dict) in sizes {
        assert_eq!(
            report(&header(&[(4, word_len)])),
            expected(1, word_len, 0, 1, 100, dict)
        );
    }
}

#[test]
fn reports_the_fields_as_the_file_holds_them() {
    // slots, step and acceleration are not checked
    let input = header(&[(11, 7), (5, u32::MAX), (6, 3), (4, 15)]);
    assert_eq!(
        report(&input),
        expected(1, 15, 7, u32::MAX, 3, "1073741824 (1073741.8k)")
    );
    // nor are the other words of the header, nor the rest of the file
    let mut input = header(&[(1, 9), (48, 9), (13, 2)]);
    input.extend(b"anything");
    assert_eq!(report(&input), expected(2, 8, 0, 1, 100, "65536 (65.5k)"));
}

#[test]
fn checks_the_number_of_sequences_against_the_size() {
    // at most one sequence for every 4 bytes of the file
    assert_eq!(
        report(&header(&[(13, 50)])),
        expected(50, 8, 0, 1, 100, "65536 (65.5k)")
    );
    let output = udbinfo(&header(&[(13, 51)]));
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "Error: invalid UDB file\n"
    );
    let mut input = header(&[(13, 51)]);
    input.extend([0; 4]);
    assert_eq!(report(&input), expected(51, 8, 0, 1, 100, "65536 (65.5k)"));
}

#[test]
fn rejects_invalid_headers() {
    for fields in [
        [(0, 0x5544_4247)],
        [(2, 64)],
        [(4, 2)],
        [(4, 16)],
        [(4, 0)],
        [(13, 0)],
        [(17, 0x0000_6161)],
        [(49, 0x5544_4246)],
    ] {
        let output = udbinfo(&header(&fields));
        assert_eq!(output.status.code(), Some(1), "{fields:?}");
        assert!(output.stdout.is_empty());
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            "Error: invalid UDB file\n"
        );
    }
}

#[test]
fn rejects_short_files() {
    for input in [&header(&[])[..199], b"", b">a\nACGT\n"] {
        let output = udbinfo(input);
        assert_eq!(output.status.code(), Some(1));
        assert!(output.stdout.is_empty());
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            "Error: invalid UDB file: shorter than its header\n"
        );
    }
}

#[test]
fn rejects_other_files() {
    // a FASTA file long enough for a header
    let fasta = format!(">a\n{}\n", "ACGT".repeat(60));
    let output = udbinfo(fasta.as_bytes());
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "Error: invalid UDB file\n"
    );
}

#[test]
fn rejects_missing_files() {
    let output = Command::new(env!("CARGO_BIN_EXE_udbinfo"))
        .args(["-i", temp_path("missing.udb").to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
}

#[test]
fn reads_and_writes_files() {
    let (input, out) = (temp_path("in.udb"), temp_path("out.txt"));
    std::fs::write(&input, header(&[(4, 5)])).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_udbinfo"))
        .args(["-i", input.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    let _ = std::fs::remove_file(&input);
    let _ = std::fs::remove_file(&out);
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(written, expected(1, 5, 0, 1, 100, "1024 (1.0k)"));
}
