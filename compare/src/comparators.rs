//! Named comparators: how a vsearch output, error or warning is compared
//! with rsearch's. Each normalisation the comparisons need lives here once.
//!
//! A file comparator takes the vsearch and rsearch files and returns `None`
//! when they agree, or a short description of the first difference. A
//! missing file reads as empty.

use std::path::Path;

fn read(path: &str) -> Vec<u8> {
    std::fs::read(path).unwrap_or_default()
}

fn lines(data: &[u8]) -> Vec<&[u8]> {
    data.split(|&b| b == b'\n').collect()
}

fn show(line: &[u8]) -> String {
    let cut = &line[..line.len().min(120)];
    format!("{:?}", String::from_utf8_lossy(cut))
}

/// The first difference between two byte strings, by line.
pub fn first_diff(v: &[u8], r: &[u8]) -> Option<String> {
    if v == r {
        return None;
    }
    let (vl, rl) = (lines(v), lines(r));
    for (k, (a, b)) in vl.iter().zip(&rl).enumerate() {
        if a != b {
            return Some(format!(
                "line {}: vsearch {} | rsearch {}",
                k + 1,
                show(a),
                show(b)
            ));
        }
    }
    Some(format!(
        "lengths differ: vsearch {} lines, rsearch {} lines",
        vl.len(),
        rl.len()
    ))
}

fn join(ls: &[&[u8]]) -> Vec<u8> {
    ls.join(&b'\n')
}

/// `--cluster_fast` writes `=` for an H record's alignment that is perfect
/// but for terminal gaps, where clust writes the CIGAR in full. An empty
/// sequence's alignment is one terminal gap.
fn perfect_but_ends(cigar: &[u8]) -> bool {
    // (\d*[ID])?\d*M(\d*[ID])? or \d*[ID]
    let s = cigar;
    let mut i = 0;
    let digits = |i: &mut usize| {
        while *i < s.len() && s[*i].is_ascii_digit() {
            *i += 1;
        }
    };
    let gap = |i: &mut usize| -> bool {
        let start = *i;
        digits(i);
        if *i < s.len() && (s[*i] == b'I' || s[*i] == b'D') {
            *i += 1;
            true
        } else {
            *i = start;
            false
        }
    };
    if gap(&mut i) && i == s.len() {
        return true;
    }
    digits(&mut i);
    if i >= s.len() || s[i] != b'M' {
        return false;
    }
    i += 1;
    gap(&mut i);
    i == s.len()
}

fn uc_fast(v: &str, r: &str) -> Option<String> {
    let (vd, rd) = (read(v), read(r));
    let (vl, rl) = (lines(&vd), lines(&rd));
    if vl.len() != rl.len() {
        return Some(format!(
            "lengths differ: vsearch {} lines, rsearch {}",
            vl.len(),
            rl.len()
        ));
    }
    for (k, (a, b)) in vl.iter().zip(&rl).enumerate() {
        let x: Vec<&[u8]> = a.split(|&c| c == b'\t').collect();
        let mut y: Vec<&[u8]> = b.split(|&c| c == b'\t').collect();
        if x.len() > 7 && y.len() > 7 && x[0] == b"H" && x[7] == b"=" && perfect_but_ends(y[7]) {
            y[7] = b"=";
        }
        if x != y {
            return Some(format!(
                "line {}: vsearch {} | rsearch {}",
                k + 1,
                show(a),
                show(b)
            ));
        }
    }
    None
}

/// BIOM without the lines naming the file, the program and the date, which
/// differ. otutab's `id` is its output file, `generated_by` names rsearch
/// and `date` is in UTC.
fn biom(data: &[u8]) -> Vec<u8> {
    let keep: Vec<&[u8]> = lines(data)
        .into_iter()
        .filter(|l| {
            !(l.starts_with(b"\t\"id\":\"")
                || contains(l, b"\"generated_by\"")
                || contains(l, b"\"date\""))
        })
        .collect();
    join(&keep)
}

fn contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

/// The report in a vsearch `--log` file, without the lines the log adds
/// before (program, command line, start time) and after (finish time,
/// elapsed time, memory).
pub fn log_report(data: &[u8]) -> Vec<u8> {
    if data.is_empty() {
        return vec![];
    }
    let ls = lines(data);
    let ls = &ls[3.min(ls.len())..];
    let Some(end) = ls.iter().rposition(|l| l.starts_with(b"Finished")) else {
        return vec![];
    };
    // the blank line before "Finished" goes too (as Python's ls[:end - 1])
    let keep = if end == 0 {
        &ls[..ls.len() - 1]
    } else {
        &ls[..end - 1]
    };
    let mut out = join(keep);
    out.push(b'\n');
    out
}

/// sed -n '/start/,/end/p': from a line matching start to the next line
/// matching end (anchored: the line starts with the pattern).
fn range(data: &[u8], start: &[u8], end: &[u8], anchored: bool) -> Vec<u8> {
    let hit = |l: &[u8], p: &[u8]| {
        if anchored {
            l.starts_with(p)
        } else {
            contains(l, p)
        }
    };
    let mut out = Vec::new();
    let mut on = false;
    for l in lines(data) {
        if !on {
            if hit(l, start) {
                on = true;
                out.extend_from_slice(l);
                out.push(b'\n');
            }
            continue;
        }
        out.extend_from_slice(l);
        out.push(b'\n');
        if hit(l, end) {
            on = false;
        }
    }
    out
}

fn grep(data: &[u8], prefix: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    for l in lines(data).into_iter().filter(|l| l.starts_with(prefix)) {
        out.extend_from_slice(l);
        out.push(b'\n');
    }
    out
}

/// `--udbstats`' report in `--log`: from the database summary (line 4) to
/// the line before the blank line and Finished.
pub fn udbstats_log(data: &[u8]) -> Vec<u8> {
    if data.is_empty() {
        return vec![];
    }
    let mut ls = lines(data);
    if ls.last() == Some(&&b""[..]) {
        ls.pop();
    }
    let mut out = Vec::new();
    if ls.len() > 7 {
        for l in &ls[3..ls.len() - 4] {
            out.extend_from_slice(l);
            out.push(b'\n');
        }
    }
    out
}

/// `--alnout` without its two opening lines, vsearch's command line and
/// version, which hitfmt does not write.
fn alnout(data: &[u8]) -> Vec<u8> {
    if data.is_empty() {
        return vec![];
    }
    let ls = lines(data);
    join(&ls[2.min(ls.len())..])
}

/// SAM without @PG lines, and the database file in @SQ's UR field.
fn sam(data: &[u8]) -> Vec<u8> {
    let mut out: Vec<Vec<u8>> = Vec::new();
    for l in lines(data) {
        if l.starts_with(b"@PG") {
            continue;
        }
        let fields: Vec<&[u8]> = l.split(|&b| b == b'\t').collect();
        let fixed: Vec<&[u8]> = fields
            .iter()
            .map(|f| {
                if f.starts_with(b"UR:file:") {
                    &b"UR:file:"[..]
                } else {
                    *f
                }
            })
            .collect();
        out.push(fixed.join(&b'\t'));
    }
    out.join(&b'\n')
}

/// Sequence lines upper-cased, FASTA or FASTQ (vsearch's `--qmask none`
/// output against `lcmask --clear`).
pub fn upper_seq(data: &[u8]) -> Vec<u8> {
    let fq = data.first() == Some(&b'@');
    let ls: Vec<Vec<u8>> = lines(data)
        .iter()
        .enumerate()
        .map(|(k, l)| {
            if if fq { k % 4 == 1 } else { !l.starts_with(b">") } {
                l.to_ascii_uppercase()
            } else {
                l.to_vec()
            }
        })
        .collect();
    ls.join(&b'\n')
}

/// FASTQ as vsearch writes it in FASTA at a width (0: unwrapped).
pub fn fastq_to_fasta(data: &[u8], width: usize) -> Vec<u8> {
    let ls = lines(data);
    let mut out = Vec::new();
    let mut k = 0;
    while k + 3 < ls.len() {
        out.push(b'>');
        out.extend_from_slice(&ls[k][1.min(ls[k].len())..]);
        out.push(b'\n');
        let s = ls[k + 1];
        if width > 0 {
            // an empty sequence: no line, as vsearch wraps
            for chunk in s.chunks(width) {
                out.extend_from_slice(chunk);
                out.push(b'\n');
            }
        } else {
            out.extend_from_slice(s);
            out.push(b'\n');
        }
        k += 4;
    }
    out
}

/// The files of `--clusters PREFIX`: the same names after the prefix, each
/// byte for byte.
fn cluster_files(v: &str, r: &str) -> Option<String> {
    let files = |prefix: &str| -> Vec<(String, String)> {
        let p = Path::new(prefix);
        let dir = p
            .parent()
            .filter(|d| !d.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let base = p
            .file_name()
            .map(|b| b.to_string_lossy().into_owned())
            .unwrap_or_default();
        let mut out: Vec<(String, String)> = std::fs::read_dir(dir)
            .map(|rd| {
                rd.filter_map(Result::ok)
                    .filter_map(|e| {
                        let name = e.file_name().to_string_lossy().into_owned();
                        name.strip_prefix(&base)
                            .map(|rest| (rest.to_string(), e.path().to_string_lossy().into_owned()))
                    })
                    .collect()
            })
            .unwrap_or_default();
        out.sort_by(|a, b| (a.0.len(), &a.0).cmp(&(b.0.len(), &b.0)));
        out
    };
    let (a, b) = (files(v), files(r));
    let keys = |x: &[(String, String)]| x.iter().map(|(k, _)| k.clone()).collect::<Vec<_>>();
    if keys(&a) != keys(&b) {
        let first5 = |k: Vec<String>| k.into_iter().take(5).collect::<Vec<_>>();
        return Some(format!(
            "files differ: vsearch {} ({:?}), rsearch {} ({:?})",
            a.len(),
            first5(keys(&a)),
            b.len(),
            first5(keys(&b))
        ));
    }
    for ((k, x), (_, y)) in a.iter().zip(&b) {
        if let Some(d) = first_diff(&read(x), &read(y)) {
            return Some(format!("cluster {k}: {d}"));
        }
    }
    None
}

// Hits that tie. vsearch sorts a query's hits (accepted before weak, by
// decreasing identity, then target) with C++'s std::sort. That leaves hits
// equal on that key in an order its library chooses. With --strand both,
// these are the two strands of one target at one identity, one match when
// the query is its own reverse complement. rsearch keeps its own order, the
// plus strand first, and keeps that strand where --maxhits keeps one. These
// comparators accept either order, and either strand for such a query.

/// The labels of the queries that are their own reverse complement, from a
/// FASTA or FASTQ file (none when it cannot be read).
fn self_complementary(path: Option<&str>) -> std::collections::HashSet<Vec<u8>> {
    let mut out = std::collections::HashSet::new();
    let Some(data) = path.map(read) else {
        return out;
    };
    let fastq = data.first() == Some(&b'@');
    let mut records: Vec<(&[u8], Vec<u8>)> = Vec::new();
    for (k, l) in lines(&data).into_iter().enumerate() {
        let l = l.strip_suffix(b"\r").unwrap_or(l);
        if fastq {
            match k % 4 {
                0 if !l.is_empty() => records.push((&l[1..], Vec::new())),
                1 => {
                    if let Some(r) = records.last_mut() {
                        r.1.extend_from_slice(l);
                    }
                }
                _ => {}
            }
        } else if let Some(label) = l.strip_prefix(b">") {
            records.push((label, Vec::new()));
        } else if let Some(r) = records.last_mut() {
            r.1.extend_from_slice(l);
        }
    }
    let complement = |c: u8| match c.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        c => c, // N, S, W
    };
    for (label, seq) in records {
        let up: Vec<u8> = seq
            .iter()
            .map(|c| match c.to_ascii_uppercase() {
                b'U' => b'T',
                c => c,
            })
            .collect();
        let rc: Vec<u8> = seq.iter().rev().map(|&c| complement(c)).collect();
        if up == rc {
            out.insert(label.to_vec());
        }
    }
    out
}

/// A hit line of a search output: its query, its tie key (none for a line
/// compared in place: a query without hits, a header), and the line with
/// its strand made neutral.
struct HitLine<'a> {
    line: &'a [u8],
    query: Vec<u8>,
    key: Option<Vec<u8>>,
    unstranded: Vec<u8>,
}

/// Reads the lines of a hit output in format `fmt`: `uc`, `blast6`, `sam`.
fn hit_lines<'a>(fmt: &str, data: &'a [u8]) -> Vec<HitLine<'a>> {
    let mut out = Vec::new();
    for line in lines(data) {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
        let tab = |fields: &[&[u8]]| fields.join(&b'\t');
        let (query, key, unstranded) = match fmt {
            // H type target length identity strand . . cigar query target
            "uc" if f.len() > 9 && f[0] == b"H" => {
                let mut u = f.clone();
                u[4] = b"*";
                (f[8].to_vec(), Some([f[1], f[3]].join(&b'\t')), tab(&u))
            }
            // query target identity ... qlo qhi: the minus strand has qlo > qhi
            "blast6" if f.len() > 7 => {
                let mut u = f.clone();
                if f[6].len() > f[7].len() || (f[6].len() == f[7].len() && f[6] > f[7]) {
                    u.swap(6, 7);
                }
                (f[0].to_vec(), Some([f[1], f[2]].join(&b'\t')), tab(&u))
            }
            // QNAME FLAG RNAME ...: the minus strand is the flag's bit 16
            "sam" if f.len() > 2 && !line.starts_with(b"@") => {
                let flag: u32 = std::str::from_utf8(f[1])
                    .ok()
                    .and_then(|t| t.parse().ok())
                    .unwrap_or(0);
                let neutral = (flag & !16).to_string();
                let mut u = f.clone();
                u[1] = neutral.as_bytes();
                let key = (flag & 4 == 0).then(|| f[2].to_vec());
                (f[0].to_vec(), key, tab(&u))
            }
            _ => (Vec::new(), None, line.to_vec()),
        };
        out.push(HitLine {
            line,
            query,
            key,
            unstranded,
        });
    }
    out
}

/// Compares two hit outputs of format `fmt`, a query's hits that tie in any
/// order. Where `--maxhits` keeps one strand of a tie, it accepts either
/// strand of a query that is its own reverse complement (`queries`, a FASTA
/// or FASTQ file).
fn hit_ties(fmt: &str, v: &str, r: &str, queries: Option<&str>) -> Option<String> {
    let (vd, rd) = (read(v), read(r));
    let (vd, rd) = if fmt == "sam" {
        (sam(&vd), sam(&rd))
    } else {
        (vd, rd)
    };
    let (vl, rl) = (hit_lines(fmt, &vd), hit_lines(fmt, &rd));
    let palindromes = self_complementary(queries);
    // blocks: runs of lines with one query and one tie key
    let blocks = |ls: &[HitLine]| -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut start = 0;
        for k in 1..=ls.len() {
            let end = k == ls.len()
                || ls[k].key.is_none()
                || ls[k].key != ls[start].key
                || ls[k].query != ls[start].query;
            if end {
                out.push((start, k));
                start = k;
            }
        }
        out
    };
    let (vb, rb) = (blocks(&vl), blocks(&rl));
    for (&(va, vz), &(ra, rz)) in vb.iter().zip(&rb) {
        let (vs, rs) = (&vl[va..vz], &rl[ra..rz]);
        let mut a: Vec<&[u8]> = vs.iter().map(|h| h.line).collect();
        let mut b: Vec<&[u8]> = rs.iter().map(|h| h.line).collect();
        a.sort();
        b.sort();
        if a == b {
            continue;
        }
        let mut au: Vec<&[u8]> = vs.iter().map(|h| h.unstranded.as_slice()).collect();
        let mut bu: Vec<&[u8]> = rs.iter().map(|h| h.unstranded.as_slice()).collect();
        au.sort();
        bu.sort();
        if au == bu && vs[0].key.is_some() && palindromes.contains(&vs[0].query) {
            continue;
        }
        return Some(format!(
            "line {}: vsearch {} | rsearch {}",
            va + 1,
            show(vs[0].line),
            show(rs[0].line)
        ));
    }
    if vb.len() != rb.len() {
        return Some(format!(
            "lengths differ: vsearch {} lines, rsearch {} lines",
            vl.len(),
            rl.len()
        ));
    }
    None
}

/// The file comparators, by name.
pub const FILES: &[&str] = &[
    "bytes",
    "uc-fast",
    "biom",
    "log-report",
    "orient-log",
    "sintax-log",
    "udbinfo-log",
    "udbstats-log",
    "alnout",
    "sam",
    "upper-seq",
    "fastq-to-fasta",
    "fastq-to-fasta0",
    "cluster-files",
    "uc-ties",
    "blast6-ties",
    "sam-ties",
];

/// Compares vsearch's file `v` with rsearch's `r` by the comparator
/// `name`: `None` when they agree.
pub fn compare_files(name: &str, v: &str, r: &str) -> Option<String> {
    // the tie comparators take the queries after a colon: uc-ties:QUERIES
    if let Some(fmt) = name.split(':').next().and_then(|n| n.strip_suffix("-ties")) {
        return hit_ties(fmt, v, r, name.split_once(':').map(|(_, q)| q));
    }
    let sides = |f: fn(&[u8]) -> Vec<u8>, both: bool| {
        let rv = read(r);
        first_diff(&f(&read(v)), &if both { f(&rv) } else { rv })
    };
    match name {
        "bytes" => first_diff(&read(v), &read(r)),
        "uc-fast" => uc_fast(v, r),
        "biom" => sides(biom, true),
        "log-report" => sides(log_report, false),
        "orient-log" => sides(
            |d| range(d, b"Forward oriented", b"Total number", true),
            false,
        ),
        "sintax-log" => sides(|d| grep(d, b"Classified"), false),
        "udbinfo-log" => sides(|d| range(d, b"Seqs", b"DBAccel", false), false),
        "udbstats-log" => sides(udbstats_log, false),
        "alnout" => sides(alnout, false),
        "sam" => sides(sam, true),
        "upper-seq" => sides(upper_seq, false),
        "fastq-to-fasta" => sides(|d| fastq_to_fasta(d, 80), false),
        "fastq-to-fasta0" => sides(|d| fastq_to_fasta(d, 0), false),
        "cluster-files" => cluster_files(v, r),
        _ => Some(format!("unknown comparator {name}")),
    }
}

// error messages: vsearch's first "Fatal error: X" against rsearch's first
// "Error: X" at the start of a line, each normalised by a kind's rules

/// The message of vsearch's first fatal error, or "".
pub fn vsearch_error(stderr: &str) -> String {
    match stderr.find("Fatal error: ") {
        Some(k) => stderr[k + 13..].lines().next().unwrap_or("").to_string(),
        None => String::new(),
    }
}

/// The message of rsearch's first error, or "".
pub fn rsearch_error(stderr: &str) -> String {
    stderr
        .lines()
        .find_map(|l| l.strip_prefix("Error: "))
        .unwrap_or("")
        .to_string()
}

/// A rewrite of a message, the effect of the regular expressions the
/// comparisons use.
enum Rule {
    /// `.*X.*`: the whole message, if it holds X.
    Holds(&'static [&'static str], &'static str),
    /// `X.*`: from X to the end.
    From(&'static str, &'static str),
    /// `.*(A|B)`: up to the last A or B.
    Upto(&'static [&'static str], &'static str),
    /// `X`: each X.
    Literal(&'static str, &'static str),
    /// `^X$`: the whole message, if it is X.
    Exact(&'static str, &'static str),
    /// `X.*(above|below).*`: from X to the end, as the last above or below
    /// after it.
    Which(&'static str),
    /// `^invalid \(([a-z]+)\) abundance annotation in read .*`.
    Joinpa,
}

fn apply(msg: &str, rules: &[Rule]) -> String {
    let mut m = msg.to_string();
    for rule in rules {
        m = match rule {
            Rule::Holds(xs, to) if xs.iter().any(|x| m.contains(x)) => to.to_string(),
            Rule::Holds(..) => m,
            Rule::From(x, to) => match m.find(x) {
                Some(k) => format!("{}{to}", &m[..k]),
                None => m,
            },
            Rule::Upto(xs, to) => match xs
                .iter()
                .filter_map(|x| m.rfind(x).map(|k| k + x.len()))
                .max()
            {
                Some(end) => format!("{to}{}", &m[end..]),
                None => m,
            },
            Rule::Literal(x, to) => m.replace(x, to),
            Rule::Exact(x, to) if m == *x => to.to_string(),
            Rule::Exact(..) => m,
            Rule::Which(x) => match m.find(x) {
                Some(k) => {
                    let rest = &m[k + x.len()..];
                    let last = ["above", "below"]
                        .iter()
                        .filter_map(|w| rest.rfind(w).map(|p| (p, *w)))
                        .max();
                    match last {
                        Some((_, w)) => format!("{}{w}", &m[..k]),
                        None => m,
                    }
                }
                None => m,
            },
            Rule::Joinpa => {
                let fixed = m.strip_prefix("invalid (").and_then(|rest| {
                    let close = rest.find(')')?;
                    let kind = &rest[..close];
                    let ok = !kind.is_empty() && kind.bytes().all(|b| b.is_ascii_lowercase());
                    (ok && rest[close..].starts_with(") abundance annotation in read ")).then(
                        || format!("Invalid ({kind}) abundance annotation in FASTA file header"),
                    )
                });
                fixed.unwrap_or(m)
            }
        };
    }
    m
}

const ABUND: Rule = Rule::Holds(&["abundance"], "abundance");

/// The two messages as an error kind compares them (`exit` compares none).
pub fn error_kinds(kind: &str, vstderr: &str, rstderr: &str) -> (String, String) {
    let (v, r) = (vsearch_error(vstderr), rsearch_error(rstderr));
    let (vr, rr): (Vec<Rule>, Vec<Rule>) = match kind {
        "invalid" => {
            // udbinfo: both say "invalid", or neither does
            let f = |m: &str| {
                if m.to_lowercase().contains("invalid") {
                    "invalid"
                } else {
                    ""
                }
            };
            return (f(&v).to_string(), f(&r).to_string());
        }
        "abundance" => (vec![ABUND], vec![ABUND]),
        "joinpa" => (vec![], vec![Rule::Joinpa]),
        "mergpa" => (
            vec![ABUND, Rule::Which("FASTQ quality value ")],
            vec![ABUND, Rule::Which("quality score ")],
        ),
        "mkudb" => (
            vec![
                ABUND,
                Rule::Holds(&["empty database", "no sequences"], "empty"),
            ],
            vec![
                ABUND,
                Rule::Holds(&["empty database", "no sequences"], "empty"),
            ],
        ),
        "mkudb-exact" => (
            vec![],
            vec![Rule::Literal(
                "no sequences: a UDB file needs at least one",
                "Cannot write a UDB file for an empty database",
            )],
        ),
        "udb" => (
            vec![
                Rule::From("Unable to read", "truncated"),
                Rule::Literal("Invalid UDB file", "invalid"),
                Rule::Literal("Incorrect UDB file size", "trailing"),
                Rule::Holds(&["too long"], "too long"),
            ],
            vec![
                Rule::Upto(&["shorter than its header", "truncated"], "truncated"),
                Rule::Upto(&["data after its end"], "trailing"),
                Rule::Exact("invalid UDB file", "invalid"),
                Rule::Holds(&["too long"], "too long"),
            ],
        ),
        _ => (vec![], vec![]), // exact
    };
    (apply(&v, &vr), apply(&r, &rr))
}

/// The error kinds a case may name.
pub const ERRORS: &[&str] = &[
    "exit",
    "exact",
    "abundance",
    "joinpa",
    "mergpa",
    "mkudb",
    "mkudb-exact",
    "udb",
    "invalid",
];

// warnings: a value extracted from each side's standard error, compared

fn count(text: &str, needle: &str) -> String {
    text.matches(needle).count().to_string()
}

fn of_sequences(text: &str, phrase: &str) -> String {
    // (\d+ of \d+) sequences PHRASE
    let tail = format!(" sequences {phrase}");
    for (k, _) in text.match_indices(&tail) {
        let head = &text[..k];
        let words: Vec<&str> = head.rsplitn(4, ' ').collect();
        if words.len() >= 3 && words[1] == "of" {
            let (a, b) = (words[2], words[0]);
            let a = a.rsplit(|c: char| !c.is_ascii_digit()).next().unwrap_or("");
            if !a.is_empty() && !b.is_empty() && b.bytes().all(|c| c.is_ascii_digit()) {
                return format!("{a} of {b}");
            }
        }
    }
    String::new()
}

fn stripped(text: &str) -> String {
    match text.find("invalid characters stripped") {
        Some(k) => text[k..].lines().next().unwrap_or("").to_string(),
        None => String::new(),
    }
}

/// The two sides of a warning, by name.
pub fn warnings(name: &str, v: &str, r: &str) -> (String, String) {
    match name {
        "stripped" => (stripped(v), stripped(r)),
        "missing-abundance" => (count(v, "Missing abundance"), count(r, "Missing abundance")),
        "nowords" => (
            of_sequences(v, "yielded no k-mer"),
            of_sequences(r, "give no word"),
        ),
        "wordlength-adjusted" => (
            count(v, "Wordlength adjusted"),
            count(r, "word length of the UDB file"),
        ),
        _ => (format!("unknown warning {name}"), String::new()),
    }
}

/// The warnings a case may name.
pub const WARNINGS: &[&str] = &[
    "stripped",
    "missing-abundance",
    "nowords",
    "wordlength-adjusted",
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigars() {
        // an empty sequence aligns as one terminal gap
        for c in ["4M", "M", "2D4M", "4MI", "3I10M2D", "339I", "D"] {
            assert!(perfect_but_ends(c.as_bytes()), "{c}");
        }
        for c in ["2MI2M", "4", "", "4M2", "2I3D", "3I3I"] {
            assert!(!perfect_but_ends(c.as_bytes()), "{c}");
        }
    }

    #[test]
    fn rules() {
        let k = |kind, v, r| error_kinds(kind, v, r);
        assert_eq!(
            k(
                "udb",
                "Fatal error: Unable to read from UDB file",
                "Error: the UDB file is truncated"
            ),
            ("truncated".into(), "truncated".into())
        );
        assert_eq!(
            k(
                "udb",
                "Fatal error: Incorrect UDB file size",
                "Error: the file has data after its end"
            ),
            ("trailing".into(), "trailing".into())
        );
        assert_eq!(
            k(
                "mergpa",
                "Fatal error: FASTQ quality value (50) above qmax (41)",
                "Error: quality score 50 is above 41"
            ),
            ("above".into(), "above".into())
        );
        assert_eq!(
            k(
                "joinpa",
                "Fatal error: Invalid (zero) abundance annotation in FASTA file header",
                "Error: invalid (zero) abundance annotation in read a;size=0"
            ),
            (
                "Invalid (zero) abundance annotation in FASTA file header".into(),
                "Invalid (zero) abundance annotation in FASTA file header".into()
            )
        );
    }

    /// Writes the vsearch, rsearch and query files of a tie test and
    /// compares them with `comparator`.
    fn ties(comparator: &str, v: &str, r: &str, queries: &str) -> Option<String> {
        let dir = std::env::temp_dir().join(format!(
            "rscompare-ties-{}-{}",
            std::process::id(),
            v.len() + 7 * r.len() + 31 * queries.len()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        let path = |n: &str| dir.join(n).to_string_lossy().into_owned();
        std::fs::write(path("v"), v).unwrap();
        std::fs::write(path("r"), r).unwrap();
        std::fs::write(path("q.fa"), queries).unwrap();
        let got = compare_files(
            &format!("{comparator}:{}", path("q.fa")),
            &path("v"),
            &path("r"),
        );
        std::fs::remove_dir_all(&dir).unwrap();
        got
    }

    const QUERIES: &str = ">q1\ncg\n>q2\nACG\n>none\nAAAA\n";

    #[test]
    fn uc_ties_in_any_order() {
        // vsearch's std::sort leaves the two strands of a target in either
        // order; rsearch writes the plus strand first
        let v = "H\t0\t2\t100.0\t-\t0\t0\t=\tq1\ts0\nH\t0\t2\t100.0\t+\t0\t0\t=\tq1\ts0\n\
                 H\t1\t2\t100.0\t+\t0\t0\t=\tq1\td1\nH\t1\t2\t100.0\t-\t0\t0\t=\tq1\td1\n\
                 N\t*\t*\t*\t.\t*\t*\t*\tnone\t*\n";
        let r = "H\t0\t2\t100.0\t+\t0\t0\t=\tq1\ts0\nH\t0\t2\t100.0\t-\t0\t0\t=\tq1\ts0\n\
                 H\t1\t2\t100.0\t+\t0\t0\t=\tq1\td1\nH\t1\t2\t100.0\t-\t0\t0\t=\tq1\td1\n\
                 N\t*\t*\t*\t.\t*\t*\t*\tnone\t*\n";
        assert_eq!(ties("uc-ties", v, r, QUERIES), None);
        // targets out of order are not a tie
        let swapped = "H\t1\t2\t100.0\t+\t0\t0\t=\tq1\td1\nH\t1\t2\t100.0\t-\t0\t0\t=\tq1\td1\n\
                       H\t0\t2\t100.0\t+\t0\t0\t=\tq1\ts0\nH\t0\t2\t100.0\t-\t0\t0\t=\tq1\ts0\n\
                       N\t*\t*\t*\t.\t*\t*\t*\tnone\t*\n";
        assert!(ties("uc-ties", v, swapped, QUERIES).is_some());
        // nor are hits of another identity
        let worse = v.replace("H\t1\t2\t100.0\t-", "H\t1\t2\t50.0\t-");
        assert!(ties("uc-ties", &worse, r, QUERIES).is_some());
        // a missing hit is a difference
        assert!(ties("uc-ties", v, &r[r.find('\n').unwrap() + 1..], QUERIES).is_some());
    }

    #[test]
    fn uc_ties_cut_by_maxhits() {
        // --maxhits 1 keeps either strand of a query that is its own
        // reverse complement. The two hits are one match
        let v = "H\t0\t2\t100.0\t-\t0\t0\t=\tq1\ts0\n";
        let r = "H\t0\t2\t100.0\t+\t0\t0\t=\tq1\ts0\n";
        assert_eq!(ties("uc-ties", v, r, QUERIES), None);
        // for any other query the strand is the hit
        let v = "H\t0\t3\t100.0\t-\t0\t0\t=\tq2\ts0\n";
        let r = "H\t0\t3\t100.0\t+\t0\t0\t=\tq2\ts0\n";
        assert!(ties("uc-ties", v, r, QUERIES).is_some());
        // and without the queries nothing is known of them
        let dir = std::env::temp_dir().join(format!("rscompare-ties-noq-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (vp, rp) = (dir.join("v"), dir.join("r"));
        std::fs::write(&vp, "H\t0\t2\t100.0\t-\t0\t0\t=\tq1\ts0\n").unwrap();
        std::fs::write(&rp, "H\t0\t2\t100.0\t+\t0\t0\t=\tq1\ts0\n").unwrap();
        let got = compare_files("uc-ties", &vp.to_string_lossy(), &rp.to_string_lossy());
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(got.is_some());
    }

    #[test]
    fn blast6_ties() {
        let v = "q1\ts0\t100.0\t2\t0\t0\t2\t1\t1\t2\t-1\t0\nq1\ts0\t100.0\t2\t0\t0\t1\t2\t1\t2\t-1\t0\n";
        let r = "q1\ts0\t100.0\t2\t0\t0\t1\t2\t1\t2\t-1\t0\nq1\ts0\t100.0\t2\t0\t0\t2\t1\t1\t2\t-1\t0\n";
        assert_eq!(ties("blast6-ties", v, r, QUERIES), None);
        // cut by --maxhits: the minus strand shows as qlo above qhi
        let v1 = "q1\ts0\t100.0\t2\t0\t0\t2\t1\t1\t2\t-1\t0\n";
        let r1 = "q1\ts0\t100.0\t2\t0\t0\t1\t2\t1\t2\t-1\t0\n";
        assert_eq!(ties("blast6-ties", v1, r1, QUERIES), None);
        let v2 = "q2\ts0\t100.0\t3\t0\t0\t3\t1\t1\t3\t-1\t0\n";
        let r2 = "q2\ts0\t100.0\t3\t0\t0\t1\t3\t1\t3\t-1\t0\n";
        assert!(ties("blast6-ties", v2, r2, QUERIES).is_some());
    }

    #[test]
    fn sam_ties() {
        let head = "@HD\tVN:1.0\tSO:unsorted\tGO:query\n@SQ\tSN:s0\tLN:2\n";
        let plus = "q1\t0\ts0\t1\t255\t2M\t*\t0\t0\tCG\t*\tAS:i:100\tMD:Z:2\n";
        let minus = "q1\t16\ts0\t1\t255\t2M\t*\t0\t0\tCG\t*\tAS:i:100\tMD:Z:2\n";
        let pg_v = "@PG\tID:vsearch\n";
        let pg_r = "@PG\tID:rsearch\n";
        assert_eq!(
            ties(
                "sam-ties",
                &format!("{head}{pg_v}{minus}{plus}"),
                &format!("{head}{pg_r}{plus}{minus}"),
                QUERIES
            ),
            None
        );
        // cut by --maxhits: the strand is the flag's bit 16
        assert_eq!(
            ties(
                "sam-ties",
                &format!("{head}{minus}"),
                &format!("{head}{plus}"),
                QUERIES
            ),
            None
        );
        let plus2 = "q2\t0\ts0\t1\t255\t3M\t*\t0\t0\tACG\t*\tAS:i:100\tMD:Z:3\n";
        let minus2 = "q2\t16\ts0\t1\t255\t3M\t*\t0\t0\tCGT\t*\tAS:i:100\tMD:Z:3\n";
        assert!(
            ties(
                "sam-ties",
                &format!("{head}{minus2}"),
                &format!("{head}{plus2}"),
                QUERIES
            )
            .is_some()
        );
        // the header still has to agree
        assert!(
            ties(
                "sam-ties",
                &format!("{head}{plus}"),
                &format!("@HD\tVN:1.0\n{plus}"),
                QUERIES
            )
            .is_some()
        );
    }

    #[test]
    fn nowords() {
        assert_eq!(
            warnings(
                "nowords",
                "WARNING: 2 of 5 sequences yielded no k-mer",
                "Warning: 2 of 5 sequences give no word"
            ),
            ("2 of 5".into(), "2 of 5".into())
        );
    }
}
