//! End-to-end tests of the sintax binary.
//!
//! Expected tables are those of `vsearch --sintax --notrunclabels --threads
//! 1` on the same files, with `--tabbedout` written to `-o`. `--seed` is
//! `--randseed`, `--cutoff` `--sintax_cutoff`, `--random` `--sintax_random`
//! and `--word-len` `--wordlength`. sintax masks the lower case of queries
//! and database alike. vsearch masks only the database's (`--dbmask dust`
//! or `soft`, without running DUST). Its expectations therefore take the
//! queries with lower case turned into N. vsearch's length filters on the
//! database (`--minseqlength 32`, `--maxseqlength 50000`) are sqfilt before
//! sintax.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

/// A path in the temporary directory, unique to this test process.
fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("rsearch-sintax-{}-{name}", std::process::id()))
}

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

    fn seq(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| b"ACGT"[(self.next() >> 62) as usize] as char)
            .collect()
    }

    /// changes one base in `one_in` at random.
    fn mutate(&mut self, seq: &str, one_in: u64) -> String {
        seq.chars()
            .map(|c| {
                if self.below(one_in) == 0 {
                    b"ACGT"[(self.next() >> 62) as usize] as char
                } else {
                    c
                }
            })
            .collect()
    }
}

/// The reverse complement of an upper-case sequence.
fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect()
}

const TAXONOMY: [Option<&str>; 6] = [
    Some(
        "d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,\
         g:Streptococcus,s:S_mutans",
    ),
    Some("d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Lactococcus"),
    Some("d:Bacteria,p:Firmicutes,c:Bacilli,o:Bacillales,f:Bacillaceae,g:Bacillus,s:B_subtilis"),
    Some("d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,g:Escherichia,s:E_coli"),
    Some("d:Archaea,p:Euryarchaeota"),
    None,
];

/// Six templates of 600 bases; a database of five copies of each, one base
/// in 25 changed, annotated with the template's taxonomy (but the last);
/// and queries of 250 bases cut from the templates, one base in 40
/// changed, one in three reverse-complemented, plus a random one, a short
/// one, a tiny one and one in lower case.
struct Data {
    templates: Vec<String>,
    db: String,
    queries: String,
}

fn data() -> Data {
    let mut rng = Rng::new(7);
    let templates: Vec<String> = (0..6).map(|_| rng.seq(600)).collect();
    let mut db = String::new();
    for i in 0..30 {
        let t = i % 6;
        let seq = rng.mutate(&templates[t], 25);
        let header = match TAXONOMY[t] {
            Some(tax) if i == 7 => format!("r{i};x=1;tax={tax}"),
            Some(tax) => format!("r{i};tax={tax};"),
            None => format!("r{i}"),
        };
        db.push_str(&format!(">{header}\n{seq}\n"));
    }
    let mut queries = String::new();
    for i in 0..16 {
        let t = i % 6;
        let start = rng.below(350) as usize;
        let mut seq = rng.mutate(&templates[t][start..start + 250], 40);
        if i % 3 == 1 {
            seq = rc(&seq);
        }
        queries.push_str(&format!(">q{i}\n{seq}\n"));
    }
    queries.push_str(&format!(">rand\n{}\n", rng.seq(250)));
    queries.push_str(&format!(">short\n{}\n", &templates[0][..35]));
    queries.push_str(">tiny\nACGT\n");
    queries.push_str(&format!(">low\n{}\n", templates[2][..200].to_lowercase()));
    Data {
        templates,
        db,
        queries,
    }
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

/// Runs sintax with `args` on `queries` against the database `db` (the
/// bytes of a FASTA or UDB file); returns the run and its report.
fn sintax(name: &str, db: &[u8], args: &[&str], queries: &str) -> (Output, String) {
    let (db_path, report) = (
        temp_path(&format!("{name}-db")),
        temp_path(&format!("{name}-report")),
    );
    std::fs::write(&db_path, db).unwrap();
    let mut all = vec![
        "--db",
        db_path.to_str().unwrap(),
        "--report",
        report.to_str().unwrap(),
    ];
    all.extend(args);
    let output = run(env!("CARGO_BIN_EXE_sintax"), &all, queries.as_bytes());
    let text = std::fs::read_to_string(&report).unwrap_or_default();
    let _ = std::fs::remove_file(&db_path);
    let _ = std::fs::remove_file(&report);
    (output, text)
}

/// Runs sintax, which must succeed, and returns its table.
fn table(name: &str, db: &[u8], args: &[&str], queries: &str) -> String {
    let (output, _) = sintax(name, db, args, queries);
    assert!(
        output.status.success(),
        "sintax failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

const TABLE: &str = r##"q0	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q1	d:Bacteria(0.77),p:Firmicutes(0.77),c:Bacilli(0.77),o:Lactobacillales(0.65),f:Streptococcaceae(0.65),g:Streptococcus(0.61),s:S_mutans(0.61)	+
q2	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q3	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q4		
q5		+
q6	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q7		
q8	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q9	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q10		
q11		+
q12	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q13		
q14	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q15	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
rand		
short		
tiny		
low		
"##;

#[test]
fn classifies_queries() {
    // the lower-case query is masked whole and not classified
    let data = data();
    let (output, report) = sintax("plus", data.db.as_bytes(), &["--seed", "1"], &data.queries);
    assert!(output.status.success());
    assert_eq!(String::from_utf8(output.stdout).unwrap(), TABLE);
    assert_eq!(report, "Classified 12 of 20 sequences (60.00%)\n");
}

#[test]
fn classifies_both_strands_with_a_cutoff() {
    let data = data();
    let args = ["--seed", "1", "--strand", "both", "--cutoff", "0.8"];
    assert_eq!(
        table("both", data.db.as_bytes(), &args, &data.queries),
        r##"q0	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+	d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus,s:S_mutans
q1	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Lactococcus(1.00)	-	d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Lactococcus
q2	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+	d:Bacteria,p:Firmicutes,c:Bacilli,o:Bacillales,f:Bacillaceae,g:Bacillus,s:B_subtilis
q3	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+	d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,g:Escherichia,s:E_coli
q4	d:Archaea(1.00),p:Euryarchaeota(1.00)	-	d:Archaea,p:Euryarchaeota
q5		+	
q6	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+	d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus,s:S_mutans
q7	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Lactococcus(1.00)	-	d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Lactococcus
q8	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+	d:Bacteria,p:Firmicutes,c:Bacilli,o:Bacillales,f:Bacillaceae,g:Bacillus,s:B_subtilis
q9	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+	d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,g:Escherichia,s:E_coli
q10	d:Archaea(1.00),p:Euryarchaeota(1.00)	-	d:Archaea,p:Euryarchaeota
q11		+	
q12	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+	d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Streptococcus,s:S_mutans
q13	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Lactococcus(1.00)	-	d:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Streptococcaceae,g:Lactococcus
q14	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+	d:Bacteria,p:Firmicutes,c:Bacilli,o:Bacillales,f:Bacillaceae,g:Bacillus,s:B_subtilis
q15	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+	d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,g:Escherichia,s:E_coli
rand	d:Bacteria(0.40),p:Firmicutes(0.37),c:Bacilli(0.37),o:Bacillales(0.21),f:Bacillaceae(0.21),g:Bacillus(0.21),s:B_subtilis(0.21)	-	
short			
tiny			
low			
"##
    );
}

#[test]
fn breaks_ties_at_random() {
    let data = data();
    let args = ["--seed", "1", "--random"];
    assert_eq!(
        table("random", data.db.as_bytes(), &args, &data.queries),
        r##"q0	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q1	d:Bacteria(0.69),p:Firmicutes(0.69),c:Bacilli(0.69),o:Lactobacillales(0.56),f:Streptococcaceae(0.56),g:Streptococcus(0.48),s:S_mutans(0.48)	+
q2	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q3	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q4		
q5		+
q6	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q7		
q8	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q9	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q10		
q11		+
q12	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q13		
q14	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q15	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
rand	d:Bacteria(0.88),p:Firmicutes(0.78),c:Bacilli(0.78),o:Lactobacillales(0.78),f:Streptococcaceae(0.78),g:Streptococcus(0.61),s:S_mutans(0.61)	+
short		
tiny		
low		
"##
    );
}

#[test]
fn takes_the_seed() {
    let data = data();
    assert_eq!(
        table("seed", data.db.as_bytes(), &["--seed", "2"], &data.queries),
        r##"q0	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q1	d:Bacteria(0.73),p:Firmicutes(0.73),c:Bacilli(0.73),o:Lactobacillales(0.59),f:Streptococcaceae(0.59),g:Streptococcus(0.47),s:S_mutans(0.47)	+
q2	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q3	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q4		
q5		+
q6	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q7		
q8	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q9	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q10		
q11		+
q12	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q13		
q14	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q15	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
rand	d:Bacteria(0.92),p:Firmicutes(0.80),c:Bacilli(0.80),o:Lactobacillales(0.80),f:Streptococcaceae(0.80),g:Streptococcus(0.72),s:S_mutans(0.72)	+
short		
tiny		
low		
"##
    );
}

#[test]
fn takes_the_word_length() {
    let data = data();
    let args = ["--seed", "1", "--word-len", "6"];
    assert_eq!(
        table("len6", data.db.as_bytes(), &args, &data.queries),
        r##"q0	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q1	d:Bacteria(0.62),p:Firmicutes(0.58),c:Bacilli(0.58),o:Lactobacillales(0.48),f:Streptococcaceae(0.48),g:Streptococcus(0.39),s:S_mutans(0.39)	+
q2	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q3	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q4	d:Bacteria(0.55),p:Firmicutes(0.47),c:Bacilli(0.47),o:Lactobacillales(0.27),f:Streptococcaceae(0.27),g:Streptococcus(0.18),s:S_mutans(0.18)	+
q5		+
q6	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q7	d:Bacteria(0.66),p:Firmicutes(0.57),c:Bacilli(0.57),o:Lactobacillales(0.33),f:Streptococcaceae(0.33),g:Streptococcus(0.19),s:S_mutans(0.19)	+
q8	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q9	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q10	d:Bacteria(0.57),p:Firmicutes(0.41),c:Bacilli(0.41),o:Lactobacillales(0.23),f:Streptococcaceae(0.23),g:Streptococcus(0.14),s:S_mutans(0.14)	+
q11		+
q12	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q13	d:Bacteria(0.45),p:Firmicutes(0.43),c:Bacilli(0.43),o:Lactobacillales(0.38),f:Streptococcaceae(0.38),g:Lactococcus(0.22)	+
q14	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q15	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
rand	d:Bacteria(0.86),p:Firmicutes(0.50),c:Bacilli(0.50),o:Lactobacillales(0.48),f:Streptococcaceae(0.48),g:Streptococcus(0.38),s:S_mutans(0.38)	+
short		
tiny		
low		
"##
    );
}

#[test]
fn skips_lowercase_database_words() {
    // the first 400 bases of the copies of the first template masked
    let data = data();
    let mut db = String::new();
    for (i, line) in data.db.lines().enumerate() {
        if i % 12 == 1 {
            db.push_str(&format!("{}{}\n", line[..400].to_lowercase(), &line[400..]));
        } else {
            db.push_str(&format!("{line}\n"));
        }
    }
    assert_eq!(
        table("lower", db.as_bytes(), &["--seed", "1"], &data.queries),
        r##"q0		
q1		
q2	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q3	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q4		
q5		+
q6	d:Bacteria(0.64),p:Firmicutes(0.46),c:Bacilli(0.46),o:Lactobacillales(0.27),f:Streptococcaceae(0.27),g:Lactococcus(0.27)	+
q7		
q8	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q9	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q10		
q11		+
q12	d:Bacteria(1.00),p:Firmicutes(0.96),c:Bacilli(0.96),o:Lactobacillales(0.88),f:Streptococcaceae(0.88),g:Streptococcus(0.85),s:S_mutans(0.85)	+
q13		
q14	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q15	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
rand		
short		
tiny		
low		
"##
    );
}

#[test]
fn reads_udb_databases() {
    let data = data();
    let udb = run(env!("CARGO_BIN_EXE_mkudb"), &[], data.db.as_bytes());
    assert!(udb.status.success());
    assert_eq!(
        table("udb", &udb.stdout, &["--seed", "1"], &data.queries),
        TABLE
    );
}

#[test]
fn parses_taxonomy() {
    // a database of one annotated template: upper-case level letters,
    // unknown ones, missing colons and empty names are skipped, a later
    // name of a level wins, and only the first tax= counts
    let data = data();
    let query = format!(">q\n{}\n", &data.templates[0][100..400]);
    let cases = [
        (
            "one;tax=D:Bac,K:Kin,z:zz,p,c:,o:Ord,o:Ord2,t:strain1",
            r##"q	d:Bac(1.00),k:Kin(1.00),o:Ord2(1.00),t:strain1(1.00)	+	d:Bac,k:Kin,o:Ord2,t:strain1
"##,
        ),
        (
            "one;tax=g:G1",
            r##"q	g:G1(1.00)	+	g:G1
"##,
        ),
        (
            "one;tax=",
            r##"q		+	
"##,
        ),
        (
            "one;tax=d:A,,p:B",
            r##"q	d:A(1.00),p:B(1.00)	+	d:A,p:B
"##,
        ),
        (
            "one;tax=d:A,p:B;tax=d:C",
            r##"q	d:A(1.00),p:B(1.00)	+	d:A,p:B
"##,
        ),
        (
            "one;tax=s:sp;x=a,b",
            r##"q	s:sp(1.00)	+	s:sp
"##,
        ),
    ];
    for (header, expected) in cases {
        let db = format!(">{header}\n{}\n", data.templates[0]);
        let args = ["--seed", "1", "--cutoff", "0.5"];
        assert_eq!(
            table("tax", db.as_bytes(), &args, &query),
            expected,
            "{header}"
        );
    }
}

#[test]
fn classifies_nothing_against_an_empty_database() {
    let data = data();
    let expected: String = data
        .queries
        .lines()
        .filter_map(|l| l.strip_prefix('>'))
        .map(|h| format!("{h}\t\t\n"))
        .collect();
    assert_eq!(
        table("empty", b"", &["--seed", "1"], &data.queries),
        expected
    );
}

#[test]
fn reports_no_queries() {
    let data = data();
    let (output, report) = sintax("none", data.db.as_bytes(), &["--seed", "1"], "");
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(report, "Classified 0 of 0 sequences\n");
}

#[test]
fn checks_abundances() {
    // the queries before an invalid abundance are written
    let data = data();
    let (output, _) = sintax(
        "zero",
        data.db.as_bytes(),
        &["--seed", "1"],
        ">a\nACGT\n>b;size=0\nACGT\n",
    );
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(String::from_utf8(output.stdout).unwrap(), "a\t\t\n");
    // the database is read whole first
    let (output, _) = sintax("zero-db", b">a;size=0\nACGT\n", &[], &data.queries);
    assert_eq!(output.status.code(), Some(1));
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("abundance"), "{stderr}");
}

#[test]
fn rejects_invalid_arguments() {
    let data = data();
    for args in [
        &["--cutoff", "1.5"][..],
        &["--cutoff", "-0.1"],
        &["--word-len", "2"],
        &["--word-len", "16"],
        &["--strand", "minus"],
    ] {
        let (output, _) = sintax("args", data.db.as_bytes(), args, &data.queries);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
    }
    // --db is required
    let output = run(env!("CARGO_BIN_EXE_sintax"), &[], b"");
    assert_eq!(output.status.code(), Some(2));
}

#[test]
fn reads_and_writes_files() {
    let data = data();
    let (db, input, out) = (
        temp_path("f-db.fa"),
        temp_path("f-in.fa"),
        temp_path("f-out"),
    );
    std::fs::write(&db, &data.db).unwrap();
    std::fs::write(&input, &data.queries).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_sintax"))
        .args(["--db", db.to_str().unwrap(), "-i", input.to_str().unwrap()])
        .args(["-o", out.to_str().unwrap(), "--seed", "1"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    let written = std::fs::read_to_string(&out).unwrap();
    for path in [&db, &input, &out] {
        let _ = std::fs::remove_file(path);
    }
    assert!(output.status.success());
    assert!(output.stdout.is_empty());
    assert_eq!(written, TABLE);
}

#[test]
fn pipeline_sintax() {
    // vsearch --sintax q --db db --randseed 1 --strand both --tabbedout t.
    // Its --minseqlength 32 drops the short database sequence, which sqfilt
    // does. It never masks queries, and lcmask --clear unmasks them
    let data = data();
    let db = format!(
        "{}>short;tax=d:Short\n{}\n",
        data.db,
        &data.templates[0][..20]
    );
    let filtered = run(
        env!("CARGO_BIN_EXE_sqfilt"),
        &["--min-len", "32", "--max-len", "50000"],
        db.as_bytes(),
    );
    assert!(filtered.status.success());
    let cleared = run(
        env!("CARGO_BIN_EXE_lcmask"),
        &["--clear"],
        data.queries.as_bytes(),
    );
    assert!(cleared.status.success());
    let args = ["--seed", "1", "--strand", "both"];
    let written = table(
        "pipe",
        &filtered.stdout,
        &args,
        std::str::from_utf8(&cleared.stdout).unwrap(),
    );
    assert_eq!(
        written,
        r##"q0	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q1	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Lactococcus(1.00)	-
q2	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q3	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q4	d:Archaea(1.00),p:Euryarchaeota(1.00)	-
q5		+
q6	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q7	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Lactococcus(1.00)	-
q8	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q9	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
q10	d:Archaea(1.00),p:Euryarchaeota(1.00)	-
q11		+
q12	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Streptococcus(1.00),s:S_mutans(1.00)	+
q13	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Lactobacillales(1.00),f:Streptococcaceae(1.00),g:Lactococcus(1.00)	-
q14	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
q15	d:Bacteria(1.00),p:Proteobacteria(1.00),c:Gammaproteobacteria(1.00),g:Escherichia(1.00),s:E_coli(1.00)	+
rand	d:Bacteria(0.40),p:Firmicutes(0.37),c:Bacilli(0.37),o:Bacillales(0.21),f:Bacillaceae(0.21),g:Bacillus(0.21),s:B_subtilis(0.21)	-
short		
tiny		
low	d:Bacteria(1.00),p:Firmicutes(1.00),c:Bacilli(1.00),o:Bacillales(1.00),f:Bacillaceae(1.00),g:Bacillus(1.00),s:B_subtilis(1.00)	+
"##
    );
}
