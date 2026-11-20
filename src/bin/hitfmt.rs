use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fmt::Write as _;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, Command, value_parser};
use md5::{Digest, Md5};

use rsearch::alignment::{Stats, parse_cigar, row};
use rsearch::cluster::{normalize_base, reverse_complement_keeping_case};
use rsearch::fastx::{Reader, Record};
use rsearch::hits::UcQueries;
use rsearch::scoring::{FOUR_BIT, N, Scoring, is_ambiguous};
use rsearch::udb::{self, Udb};

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// vsearch's user fields, in its order.
const FIELDS: [&str; 53] = [
    "query", "target", "evalue", "id", "pctpv", "pctgaps", "pairs", "gaps", "qlo", "qhi", "tlo",
    "thi", "pv", "ql", "tl", "qs", "ts", "alnlen", "opens", "exts", "raw", "bits", "aln", "caln",
    "qstrand", "tstrand", "qrow", "trow", "qframe", "tframe", "mism", "ids", "qcov", "tcov", "id0",
    "id1", "id2", "id3", "id4", "qilo", "qihi", "tilo", "tihi", "diffs", "mid", "qseq", "tseq",
    "qrowdots", "trowdots", "qlor", "qhir", "tlor", "thir",
];

/// parses user fields separated by `+`.
fn user_fields(text: &str) -> Result<Vec<usize>, String> {
    text.split('+')
        .map(|name| {
            FIELDS
                .iter()
                .position(|&f| f == name)
                .ok_or_else(|| format!("unknown user field '{name}'"))
        })
        .collect()
}

/// The output format.
enum Format {
    Alignments,
    Blast6,
    User(Vec<usize>),
    Sam,
    FastaPairs,
    QuerySegments,
    TargetSegments,
}

/// A hit, rebuilt.
struct Hit {
    target: usize,
    minus: bool,
    cigar: String,
    ops: Vec<(u8, usize)>,
    stats: Stats,
}

/// A query: its header, its sequence and its reverse complement.
struct Query {
    header: String,
    seq: Vec<u8>,
    rc: Vec<u8>,
}

impl Query {
    fn strand(&self, minus: bool) -> &[u8] {
        if minus { &self.rc } else { &self.seq }
    }
}

struct Formatter<'a> {
    format: Format,
    db: &'a [(String, Vec<u8>)],
    scoring: Scoring,
    iddef: usize,
    rowlen: usize,
    no_hits: bool,
    /// vsearch's alignment symbol of each pair of 4-bit codes.
    symbols: [[u8; 16]; 16],
}

/// formats a share with one decimal, as `%.1f`.
fn one(x: f64) -> String {
    format!("{x:.1}")
}

impl Formatter<'_> {
    fn id(&self, hit: &Hit) -> f64 {
        hit.stats.id[self.iddef]
    }

    /// a row of the hit's alignment, without its terminal gaps.
    fn window(&self, hit: &Hit, row: &[u8]) -> Vec<u8> {
        let s = &hit.stats;
        let from = (s.trim_q_left + s.trim_t_left) as usize;
        row[from..from + s.internal_length as usize].to_vec()
    }

    fn rows(&self, query: &Query, hit: &Hit) -> (Vec<u8>, Vec<u8>) {
        let target = &self.db[hit.target].1;
        (
            row(query.strand(hit.minus), &hit.ops, true),
            row(target, &hit.ops, false),
        )
    }

    fn write(&self, out: &mut Vec<u8>, query: &Query, hits: &[Hit]) {
        match &self.format {
            Format::Alignments => self.alignments(out, query, hits),
            Format::Blast6 => {
                if hits.is_empty() && self.no_hits {
                    out.extend_from_slice(query.header.as_bytes());
                    out.extend_from_slice(b"\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n");
                }
                for hit in hits {
                    let s = &hit.stats;
                    let (qlen, tlen) = (query.seq.len(), self.db[hit.target].1.len());
                    let (qstart, qend) = if hit.minus { (qlen, 1) } else { (1, qlen) };
                    let _ = writeln!(
                        Bytes(out),
                        "{}\t{}\t{}\t{}\t{}\t{}\t{qstart}\t{qend}\t1\t{tlen}\t-1\t0",
                        query.header,
                        self.db[hit.target].0,
                        one(self.id(hit)),
                        s.internal_length,
                        s.mismatches,
                        s.internal_gaps
                    );
                }
            }
            Format::User(fields) => {
                if hits.is_empty() && self.no_hits {
                    self.user_line(out, fields, query, None);
                }
                for hit in hits {
                    self.user_line(out, fields, query, Some(hit));
                }
            }
            Format::Sam => self.sam(out, query, hits),
            Format::FastaPairs => {
                for hit in hits {
                    let (qrow, trow) = self.rows(query, hit);
                    fasta(out, &query.header, &self.window(hit, &qrow));
                    fasta(out, &self.db[hit.target].0, &self.window(hit, &trow));
                    out.push(b'\n');
                }
            }
            Format::QuerySegments => {
                for hit in hits {
                    let seq = query.strand(hit.minus);
                    let s = &hit.stats;
                    let end = seq.len() as i64 - s.trim_q_right;
                    fasta(
                        out,
                        &query.header,
                        &seq[s.trim_q_left as usize..end as usize],
                    );
                }
            }
            Format::TargetSegments => {
                for hit in hits {
                    let (header, seq) = &self.db[hit.target];
                    let s = &hit.stats;
                    let end = seq.len() as i64 - s.trim_t_right;
                    fasta(out, header, &seq[s.trim_t_left as usize..end as usize]);
                }
            }
        }
    }

    fn alignments(&self, out: &mut Vec<u8>, query: &Query, hits: &[Hit]) {
        if hits.is_empty() {
            if self.no_hits {
                let _ = write!(Bytes(out), "\nQuery >{}\nNo hits\n", query.header);
            }
            return;
        }
        let qlen = query.seq.len();
        let _ = write!(
            Bytes(out),
            "\nQuery >{}\n %Id   TLen  Target\n",
            query.header
        );
        for hit in hits {
            let (header, seq) = &self.db[hit.target];
            let _ = writeln!(
                Bytes(out),
                "{:3.0}% {:6}  {header}",
                self.id(hit),
                seq.len()
            );
        }
        for hit in hits {
            let (header, target) = &self.db[hit.target];
            let tlen = target.len();
            let width = qlen.to_string().len().max(tlen.to_string().len());
            let _ = write!(
                Bytes(out),
                "\n Query {qlen:width$}nt >{}\nTarget {tlen:width$}nt >{header}\n",
                query.header
            );
            let rowlen = if self.rowlen == 0 {
                qlen + tlen
            } else {
                self.rowlen
            };
            self.show(out, query, hit, target, width, rowlen);
            let s = &hit.stats;
            let gaps = if s.internal_length > 0 {
                100.0 * s.internal_indels as f64 / s.internal_length as f64
            } else {
                0.0
            };
            let _ = write!(
                Bytes(out),
                "\n{} cols, {} ids ({:3.1}%), {} gaps ({gaps:3.1}%)\n",
                s.internal_length,
                s.matches,
                self.id(hit),
                s.internal_indels
            );
        }
    }

    /// vsearch's `align_show`: blocks of `rowlen` columns of the alignment
    /// without its terminal gaps, the query on the strand it was aligned,
    /// numbered on the plus strand.
    fn show(
        &self,
        out: &mut Vec<u8>,
        query: &Query,
        hit: &Hit,
        target: &[u8],
        width: usize,
        rowlen: usize,
    ) {
        let s = &hit.stats;
        let qlen = query.seq.len() as i64;
        let qseq = query.strand(hit.minus);
        let cigar = &hit.cigar[s.trim_aln_left..hit.cigar.len() - s.trim_aln_right];
        // the three rows, whole
        let (mut qrow, mut trow) = (Vec::new(), Vec::new());
        let (mut q, mut t) = (s.trim_q_left as usize, s.trim_t_left as usize);
        for (op, run) in parse_cigar(cigar).unwrap() {
            match op {
                b'M' => {
                    qrow.extend_from_slice(&qseq[q..q + run]);
                    trow.extend_from_slice(&target[t..t + run]);
                    q += run;
                    t += run;
                }
                b'D' => {
                    qrow.extend_from_slice(&qseq[q..q + run]);
                    trow.extend(std::iter::repeat_n(b'-', run));
                    q += run;
                }
                _ => {
                    qrow.extend(std::iter::repeat_n(b'-', run));
                    trow.extend_from_slice(&target[t..t + run]);
                    t += run;
                }
            }
        }
        let symbols: Vec<u8> = qrow
            .iter()
            .zip(&trow)
            .map(|(&a, &b)| {
                if a == b'-' || b == b'-' {
                    b' '
                } else {
                    self.symbols[usize::from(FOUR_BIT[usize::from(a)])]
                        [usize::from(FOUR_BIT[usize::from(b)])]
                }
            })
            .collect();
        // vsearch's positions: the query's index on the plus strand, moving
        // down on the minus strand
        let (mut qpos, delta) = if hit.minus {
            (qlen - 1 - s.trim_q_left, -1)
        } else {
            (s.trim_q_left, 1)
        };
        let mut tpos = s.trim_t_left;
        let rowlen = rowlen.max(1);
        for start in (0..qrow.len()).step_by(rowlen) {
            let end = (start + rowlen).min(qrow.len());
            let (qstart, tstart) = (qpos, tpos);
            let residues = |row: &[u8]| row.iter().filter(|&&c| c != b'-').count() as i64;
            qpos += delta * residues(&qrow[start..end]);
            tpos += residues(&trow[start..end]);
            let query_start = (qstart + 1).min(qlen);
            let query_end = if hit.minus { qpos + 2 } else { qpos };
            let target_start = (tstart + 1).min(target.len() as i64);
            let _ = write!(
                Bytes(out),
                "\nQry {query_start:width$} {} ",
                if hit.minus { '-' } else { '+' }
            );
            out.extend_from_slice(&qrow[start..end]);
            let _ = write!(Bytes(out), " {query_end}\n    {:width$}   ", "");
            out.extend_from_slice(&symbols[start..end]);
            let _ = write!(Bytes(out), "\nTgt {target_start:width$} + ");
            out.extend_from_slice(&trow[start..end]);
            let _ = writeln!(Bytes(out), " {tpos}");
        }
    }

    fn sam(&self, out: &mut Vec<u8>, query: &Query, hits: &[Hit]) {
        if hits.is_empty() {
            if self.no_hits {
                let _ = writeln!(
                    Bytes(out),
                    "{}\t4\t*\t0\t255\t*\t*\t0\t0\t{}\t*",
                    query.header,
                    String::from_utf8_lossy(&query.seq)
                );
            }
            return;
        }
        for (t, hit) in hits.iter().enumerate() {
            let s = &hit.stats;
            let qseq = query.strand(hit.minus);
            let (header, target) = &self.db[hit.target];
            let (cigar, md) = sam_strings(&hit.ops, qseq, target);
            let flag = (0x10 * u32::from(hit.minus)) | if t > 0 { 0x100 } else { 0 };
            let _ = writeln!(
                Bytes(out),
                "{}\t{flag}\t{header}\t1\t255\t{cigar}\t*\t0\t0\t{}\t*\tAS:i:{:.0}\tXN:i:0\tXM:i:{}\tXO:i:{}\tXG:i:{}\tNM:i:{}\tMD:Z:{md}\tYT:Z:UU",
                query.header,
                String::from_utf8_lossy(qseq),
                self.id(hit),
                s.mismatches,
                s.internal_gaps,
                s.internal_indels,
                s.mismatches + s.internal_indels
            );
        }
    }

    fn user_line(&self, out: &mut Vec<u8>, fields: &[usize], query: &Query, hit: Option<&Hit>) {
        for (i, &field) in fields.iter().enumerate() {
            if i > 0 {
                out.push(b'\t');
            }
            self.user_field(out, FIELDS[field], query, hit);
        }
        out.push(b'\n');
    }

    fn user_field(&self, out: &mut Vec<u8>, field: &str, query: &Query, hit: Option<&Hit>) {
        let qlen = query.seq.len() as i64;
        let Some(hit) = hit else {
            let text = match field {
                "query" => query.header.clone(),
                "target" => "*".into(),
                "evalue" => "-1".into(),
                "id" | "pctpv" | "pctgaps" | "qcov" | "tcov" | "id0" | "id1" | "id2" | "id3"
                | "id4" | "mid" => "0.0".into(),
                "ql" | "qs" => qlen.to_string(),
                "qframe" | "tframe" => "+0".into(),
                "qseq" => String::from_utf8_lossy(&query.seq).into_owned(),
                "aln" | "caln" | "qstrand" | "tstrand" | "qrow" | "trow" | "tseq" | "qrowdots"
                | "trowdots" => String::new(),
                _ => "0".into(),
            };
            out.extend_from_slice(text.as_bytes());
            return;
        };
        let s = &hit.stats;
        let (theader, tseq) = &self.db[hit.target];
        let tlen = tseq.len() as i64;
        let pct = |x: i64, of: i64| {
            if of > 0 {
                100.0 * x as f64 / of as f64
            } else {
                0.0
            }
        };
        let text: String = match field {
            "query" => query.header.clone(),
            "target" => theader.clone(),
            "evalue" => "-1".into(),
            "id" => one(self.id(hit)),
            "pctpv" => one(pct(s.matches, s.internal_length)),
            "pctgaps" => one(pct(s.internal_indels, s.internal_length)),
            "pairs" => (s.matches + s.mismatches).to_string(),
            "gaps" => s.internal_indels.to_string(),
            "qlo" => (if hit.minus { qlen } else { 1 }).to_string(),
            "qhi" => (if hit.minus { 1 } else { qlen }).to_string(),
            "tlo" => "1".into(),
            "thi" | "tl" | "ts" => tlen.to_string(),
            "pv" | "ids" => s.matches.to_string(),
            "ql" | "qs" => qlen.to_string(),
            "alnlen" => s.internal_length.to_string(),
            "opens" => s.internal_gaps.to_string(),
            "exts" => (s.internal_indels - s.internal_gaps).to_string(),
            "raw" => s.score.to_string(),
            "bits" => "0".into(),
            "aln" => hit
                .ops
                .iter()
                .flat_map(|&(op, run)| std::iter::repeat_n(op as char, run))
                .collect(),
            "caln" => hit.cigar.clone(),
            "qstrand" => (if hit.minus { "-" } else { "+" }).into(),
            "tstrand" => "+".into(),
            "qrow" => {
                let (qrow, _) = self.rows(query, hit);
                String::from_utf8_lossy(&self.window(hit, &qrow)).into_owned()
            }
            "trow" => {
                let (_, trow) = self.rows(query, hit);
                String::from_utf8_lossy(&self.window(hit, &trow)).into_owned()
            }
            "qframe" | "tframe" => "+0".into(),
            "mism" => s.mismatches.to_string(),
            "qcov" => one(100.0 * (s.matches + s.mismatches) as f64 / qlen as f64),
            "tcov" => one(100.0 * (s.matches + s.mismatches) as f64 / tlen as f64),
            "id0" => one(s.id[0]),
            "id1" => one(s.id[1]),
            "id2" => one(s.id[2]),
            "id3" => one(s.id[3]),
            "id4" => one(s.id[4]),
            "qilo" => (s.trim_q_left + 1).to_string(),
            "qihi" => (qlen - s.trim_q_right).to_string(),
            "tilo" => (s.trim_t_left + 1).to_string(),
            "tihi" => (tlen - s.trim_t_right).to_string(),
            "diffs" => (s.mismatches + s.internal_indels).to_string(),
            "mid" => one(pct(s.matches, s.matches + s.mismatches)),
            "qseq" => String::from_utf8_lossy(query.strand(hit.minus)).into_owned(),
            "tseq" => String::from_utf8_lossy(tseq).into_owned(),
            "qrowdots" | "trowdots" => {
                let (qrow, trow) = self.rows(query, hit);
                let (qrow, trow) = (self.window(hit, &qrow), self.window(hit, &trow));
                let (own, other) = if field == "qrowdots" {
                    (qrow, trow)
                } else {
                    (trow, qrow)
                };
                own.iter()
                    .zip(&other)
                    .map(|(&a, &b)| {
                        let same = a != b'-'
                            && b != b'-'
                            && FOUR_BIT[usize::from(a)] == FOUR_BIT[usize::from(b)];
                        if same { '.' } else { a as char }
                    })
                    .collect()
            }
            "qlor" => (if hit.minus { qlen - 1 } else { 0 }).to_string(),
            "qhir" => (if hit.minus { 0 } else { qlen - 1 }).to_string(),
            "tlor" => "0".into(),
            "thir" => (tlen - 1).to_string(),
            _ => unreachable!("a field of FIELDS"),
        };
        out.extend_from_slice(text.as_bytes());
    }
}

/// vsearch's alignment symbol of each pair of 4-bit codes: `|` for a base
/// against itself, `+` for codes sharing a base, a space otherwise and,
/// with `--n-mismatch`, against N.
fn symbols(n_mismatch: bool) -> [[u8; 16]; 16] {
    let mut table = [[b' '; 16]; 16];
    for a in 0..16u8 {
        for b in 0..16u8 {
            table[usize::from(a)][usize::from(b)] = if n_mismatch && (a == N || b == N) {
                b' '
            } else if a == b && !is_ambiguous(a) {
                b'|'
            } else if a & b != 0 {
                b'+'
            } else {
                b' '
            };
        }
    }
    table
}

/// writes a FASTA record, unwrapped.
fn fasta(out: &mut Vec<u8>, header: &str, seq: &[u8]) {
    out.push(b'>');
    out.extend_from_slice(header.as_bytes());
    out.push(b'\n');
    out.extend_from_slice(seq);
    out.push(b'\n');
}

/// vsearch's SAM CIGAR (with I and D swapped, every run with its length)
/// and MD string of an alignment.
fn sam_strings(ops: &[(u8, usize)], query: &[u8], target: &[u8]) -> (String, String) {
    let (mut cigar, mut md) = (String::new(), String::new());
    let (mut q, mut t) = (0, 0);
    let mut matched = 0;
    let mut ends_with_number = false;
    for &(op, run) in ops {
        match op {
            b'M' => {
                let _ = write!(cigar, "{run}M");
                for _ in 0..run {
                    if FOUR_BIT[usize::from(query[q])] == FOUR_BIT[usize::from(target[t])] {
                        matched += 1;
                    } else {
                        if !ends_with_number {
                            let _ = write!(md, "{matched}");
                            matched = 0;
                        }
                        md.push(target[t].to_ascii_uppercase() as char);
                        ends_with_number = false;
                    }
                    q += 1;
                    t += 1;
                }
            }
            b'D' => {
                let _ = write!(cigar, "{run}I");
                q += run;
            }
            _ => {
                let _ = write!(cigar, "{run}D");
                if !ends_with_number {
                    let _ = write!(md, "{matched}");
                    matched = 0;
                }
                md.push('^');
                for _ in 0..run {
                    md.push(target[t].to_ascii_uppercase() as char);
                    t += 1;
                }
                ends_with_number = false;
            }
        }
    }
    if !ends_with_number {
        let _ = write!(md, "{matched}");
    }
    (cigar, md)
}

/// A `fmt::Write` over a byte buffer.
struct Bytes<'a>(&'a mut Vec<u8>);

impl std::fmt::Write for Bytes<'_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        self.0.extend_from_slice(s.as_bytes());
        Ok(())
    }
}

/// reads the database: a FASTA, FASTQ or UDB file.
fn database(path: &str) -> Vec<(String, Vec<u8>)> {
    let file = File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}")));
    let size = file
        .metadata()
        .ok()
        .filter(|m| m.is_file())
        .map(|m| m.len());
    let mut fin = BufReader::new(file);
    let start = fin
        .fill_buf()
        .unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
    if start.starts_with(udb::MAGIC) {
        let ignore = |_: &[u32], _: u64, _: &[u32]| {};
        let read = match size {
            Some(size) => Udb::read(&mut fin, size, ignore),
            None => {
                let mut bytes = Vec::new();
                fin.read_to_end(&mut bytes)
                    .unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
                Udb::read(&mut &bytes[..], bytes.len() as u64, ignore)
            }
        };
        let udb = read.unwrap_or_else(|e| fail(&e.to_string()));
        let mut at = 0;
        return udb
            .header_starts
            .windows(2)
            .zip(&udb.lengths)
            .map(|(pair, &length)| {
                let header = &udb.headers[pair[0] as usize..pair[1] as usize - 1];
                let seq = udb.sequences[at..at + length as usize].to_vec();
                at += length as usize;
                (String::from_utf8_lossy(header).into_owned(), seq)
            })
            .collect();
    }
    let mut reader = Reader::new(fin);
    let mut record = Record::new();
    let mut db = Vec::new();
    loop {
        reader
            .read(&mut record)
            .unwrap_or_else(|e| fail(&e.to_string()));
        if record.is_empty() {
            return db;
        }
        db.push((record.header().to_string(), record.seq().to_vec()));
    }
}

fn main() {
    let flag = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .action(ArgAction::SetTrue)
            .help(help)
    };
    let gaps = |text: &str, open: bool| -> Result<String, String> {
        Scoring::default()
            .set_gaps(text, open)
            .map(|()| text.to_string())
    };
    let mut command = Command::new("hitfmt")
        .about("Write hits as alignments, BLAST tabular, user fields, SAM, FASTA pairs or aligned segments.")
        .long_about(
            "Write hits as alignments, BLAST tabular, user fields, SAM, FASTA pairs or aligned segments.\n\n\
             Rebuilds vsearch's hit outputs from UC hits, the queries and the database \
             they were found in: the queries are matched to the hits in order, and \
             targets by their number. Give the scores and the identity definition the \
             search used.",
        )
        .arg(Arg::new("hits").value_name("HITS").required(true).help("UC file of hits"))
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .required(true)
                .help("FASTA or FASTQ file of the queries"),
        )
        .arg(
            Arg::new("db")
                .long("db")
                .value_name("FILE")
                .required(true)
                .help("Database: FASTA, FASTQ or UDB file"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file, default to stdout"),
        )
        .arg(flag("alignments", "Write alignments, as vsearch's --alnout"))
        .arg(flag("blast6", "Write BLAST tabular, as vsearch's --blast6out"))
        .arg(
            Arg::new("userout")
                .long("userout")
                .value_name("FIELDS")
                .value_parser(|s: &str| user_fields(s).map(|_| s.to_string()))
                .help("Write these user fields, separated by +, as vsearch's --userout"),
        )
        .arg(flag("sam", "Write SAM, as vsearch's --samout"))
        .arg(flag("fasta-pairs", "Write the aligned rows of each hit, as vsearch's --fastapairs"))
        .arg(flag("qsegs", "Write the aligned segments of the queries, as vsearch's --qsegout"))
        .arg(flag("tsegs", "Write the aligned segments of the targets, as vsearch's --tsegout"))
        .group(
            ArgGroup::new("format")
                .args([
                    "alignments",
                    "blast6",
                    "userout",
                    "sam",
                    "fasta-pairs",
                    "qsegs",
                    "tsegs",
                ])
                .required(true),
        )
        .arg(
            Arg::new("sam-header")
                .long("sam-header")
                .action(ArgAction::SetTrue)
                .help("Begin SAM with a header"),
        )
        .arg(flag(
            "output-no-hits",
            "Also write queries without hits, as vsearch's --output_no_hits",
        ))
        .arg(
            Arg::new("rowlen")
                .long("rowlen")
                .value_name("LEN")
                .value_parser(value_parser!(usize))
                .default_value("64")
                .help("Columns of an alignment row, 0 for one row"),
        )
        .arg(
            Arg::new("iddef")
                .long("iddef")
                .value_name("DEF")
                .value_parser(value_parser!(u8).range(0..=4))
                .default_value("2")
                .help("Identity definition, 0 to 4, as the search's"),
        )
        .arg(flag("n-mismatch", "Count N as a mismatch, as the search did"))
        .arg(
            Arg::new("match")
                .long("match")
                .value_name("SCORE")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i64).range(-32767..=32767))
                .default_value("2")
                .help("Score of a match, as the search's"),
        )
        .arg(
            Arg::new("mismatch")
                .long("mismatch")
                .value_name("SCORE")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i64).range(-32767..=32767))
                .default_value("-4")
                .help("Score of a mismatch, as the search's"),
        )
        .arg(
            Arg::new("gapopen")
                .long("gapopen")
                .value_name("PENALTIES")
                .value_parser(move |s: &str| gaps(s, true))
                .help("Gap opening penalties, as the search's, default to 20I/2E"),
        )
        .arg(
            Arg::new("gapext")
                .long("gapext")
                .value_name("PENALTIES")
                .value_parser(move |s: &str| gaps(s, false))
                .help("Gap extension penalties, as the search's, default to 2I/1E"),
        );
    let matches = command.get_matches_mut();
    if matches.get_flag("sam-header") && !matches.get_flag("sam") {
        command
            .error(
                clap::error::ErrorKind::MissingRequiredArgument,
                "--sam-header needs --sam",
            )
            .exit();
    }

    let n_mismatch = matches.get_flag("n-mismatch");
    let mut scoring = Scoring {
        match_score: *matches.get_one::<i64>("match").unwrap(),
        mismatch: *matches.get_one::<i64>("mismatch").unwrap(),
        n_mismatch,
        ..Scoring::default()
    };
    if let Some(text) = matches.get_one::<String>("gapopen") {
        scoring.set_gaps(text, true).unwrap();
    }
    if let Some(text) = matches.get_one::<String>("gapext") {
        scoring.set_gaps(text, false).unwrap();
    }
    let format = if matches.get_flag("alignments") {
        Format::Alignments
    } else if matches.get_flag("blast6") {
        Format::Blast6
    } else if let Some(fields) = matches.get_one::<String>("userout") {
        Format::User(user_fields(fields).unwrap())
    } else if matches.get_flag("sam") {
        Format::Sam
    } else if matches.get_flag("fasta-pairs") {
        Format::FastaPairs
    } else if matches.get_flag("qsegs") {
        Format::QuerySegments
    } else {
        Format::TargetSegments
    };

    let db_path = matches.get_one::<String>("db").unwrap();
    let db = database(db_path);
    let hits_path = matches.get_one::<String>("hits").unwrap();
    let hits_file =
        File::open(hits_path).unwrap_or_else(|e| fail(&format!("cannot open {hits_path}: {e}")));
    let queries_path = matches.get_one::<String>("in").unwrap();
    let queries_file = File::open(queries_path)
        .unwrap_or_else(|e| fail(&format!("cannot open {queries_path}: {e}")));
    let mut queries = Reader::new(BufReader::new(queries_file));
    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let formatter = Formatter {
        format,
        db: &db,
        scoring,
        iddef: usize::from(*matches.get_one::<u8>("iddef").unwrap()),
        rowlen: *matches.get_one::<usize>("rowlen").unwrap(),
        no_hits: matches.get_flag("output-no-hits"),
        symbols: symbols(n_mismatch),
    };
    let table = formatter.scoring.substitutions();
    let mut buffer = Vec::new();
    if matches.get_flag("sam-header") {
        buffer.extend_from_slice(b"@HD\tVN:1.0\tSO:unsorted\tGO:query\n");
        for (header, seq) in &db {
            let mut hasher = Md5::new();
            hasher.update(seq.iter().map(|&c| normalize_base(c)).collect::<Vec<u8>>());
            let digest: String = hasher
                .finalize()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            let _ = writeln!(
                Bytes(&mut buffer),
                "@SQ\tSN:{header}\tLN:{}\tM5:{digest}\tUR:file:{db_path}",
                seq.len()
            );
        }
        let command: Vec<String> = std::env::args().collect();
        let _ = writeln!(
            Bytes(&mut buffer),
            "@PG\tID:rsearch\tVN:{}\tCL:{}",
            env!("CARGO_PKG_VERSION"),
            command.join(" ")
        );
    }

    let mut uc = UcQueries::new(BufReader::new(hits_file));
    let mut record = Record::new();
    loop {
        let target = |t: u64| db.get(t as usize).map(|(header, _)| header.as_str());
        let (label, group) = match uc.next_query(target) {
            None => break,
            Some(result) => result.unwrap_or_else(|e| fail(&e)),
        };
        // the query, skipping those without records
        loop {
            queries
                .read(&mut record)
                .unwrap_or_else(|e| fail(&e.to_string()));
            if record.is_empty() {
                fail(&format!(
                    "query {label} of the hits not found in {queries_path}"
                ));
            }
            if record.header() == label {
                break;
            }
        }
        let mut query = Query {
            header: label.clone(),
            seq: record.seq().to_vec(),
            rc: Vec::new(),
        };
        reverse_complement_keeping_case(&mut query.rc, &query.seq);
        let hits: Vec<Hit> = group
            .into_iter()
            .map(|h| {
                if h.query_len != query.seq.len() {
                    fail(&format!("hit of query {label} of another length"));
                }
                let cigar = if h.cigar == "=" {
                    format!("{}M", query.seq.len())
                } else {
                    h.cigar
                };
                let target = h.target as usize;
                let stats = Stats::new(
                    &cigar,
                    query.strand(h.minus),
                    &db[target].1,
                    &formatter.scoring,
                    &table,
                )
                .unwrap_or_else(|e| fail(&format!("{e}: hit of query {label}")));
                Hit {
                    target,
                    minus: h.minus,
                    ops: parse_cigar(&cigar).unwrap(),
                    cigar,
                    stats,
                }
            })
            .collect();
        formatter.write(&mut buffer, &query, &hits);
        out.write_all(&buffer)
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
        buffer.clear();
    }
    out.write_all(&buffer)
        .and_then(|()| out.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
