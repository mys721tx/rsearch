//! Generators for the preprocessing tools' cases: the edge cases of each
//! tool's input (empty and short reads, lower case, IUPAC codes, labels
//! with descriptions and abundances, encodings, malformed abundances).

use std::io;
use std::path::Path;

use super::{Text, fasta, fastq, rc, w, wb};
use crate::rng::Rng;
use crate::translate::py_float;

fn size(r: &mut Rng, a: i64, b: i64) -> String {
    format!(";size={}", r.randint(a, b))
}

pub fn derep(t: &Path, r: &mut Rng) -> io::Result<()> {
    let base: Vec<String> = (0..12)
        .map(|_| {
            let n = r.randint(8, 14) as usize;
            r.chars("ACGT", n)
        })
        .collect();
    let (mut fq, mut fa, mut pre) = (Text::new(), Text::new(), Text::new());
    for i in 0..300 {
        let mut s = r.choice(&base);
        let k = r.random();
        if k < 0.15 {
            s = rc(&s);
        } else if k < 0.25 {
            s = s.to_lowercase();
        }
        if r.chance(0.03) {
            s.clear();
        }
        let mut lab = format!("r{}", r.randint(0, 40));
        if r.chance(0.5) {
            lab += &size(r, 1, 9);
        }
        if r.chance(0.2) {
            lab += r.choice(&[" d1", " d2", "\td3", ";x=1 d"]);
        }
        let q: String = s
            .bytes()
            .map(|_| {
                (33 + r.choice(&[0, 1, 2, 3, 5, 10, 17, 20, 25, 30, 35, 38, 40, 41])) as u8 as char
            })
            .collect();
        fq.push(&fastq(&lab, &s, &q));
        fa.push(&fasta(&lab, &s));
        let p = r.choice(&base);
        let lo = if r.chance(0.03) { 0 } else { 3 };
        let p = &p[..r.randint(lo, p.len() as i64) as usize];
        let size_ = if r.chance(0.5) {
            size(r, 1, 4)
        } else {
            String::new()
        };
        pre.push(&format!(">p{}{size_}\n{p}\n", i % 37));
    }
    w(t, "q.fq", &fq.0)?;
    w(t, "a.fa", &fa.0)?;
    w(t, "p.fa", &pre.0)?;
    // the same reads at Phred+64
    let q64: String =
        fq.0.split_inclusive('\n')
            .enumerate()
            .map(|(i, l)| {
                if i % 4 == 3 {
                    let body: String = l
                        .trim_end_matches('\n')
                        .bytes()
                        .map(|c| (c + 31) as char)
                        .collect();
                    body + "\n"
                } else {
                    l.to_string()
                }
            })
            .collect();
    w(t, "q64.fq", &q64)
}

fn conv_label(r: &mut Rng, i: usize) -> String {
    let mut l = format!("r{i}");
    if r.chance(0.4) {
        l += &size(r, 1, 9);
    }
    if r.chance(0.3) {
        l += ";ee=0.25";
    }
    if r.chance(0.3) {
        l += ";length=7";
    }
    if r.chance(0.3) {
        l += r.choice(&[" desc", " a b", "\tt"]);
    }
    l
}

/// Abundances that are invalid, added to three random headers (every
/// `step`th line).
fn spoil(r: &mut Rng, text: &str, step: usize) -> String {
    let mut lines: Vec<String> = text.split('\n').map(str::to_string).collect();
    let records = lines.len() / step;
    if records > 1 {
        for k in r.sample(records - 1, 3) {
            let k = k + 1;
            lines[step * k] += r.choice(&[";size=0", ";size=9223372036854775808"]);
        }
    }
    lines.join("\n")
}

pub fn conv(t: &Path, r: &mut Rng) -> io::Result<()> {
    let (mut q33, mut q64, mut sol, mut fa) = (Text::new(), Text::new(), Text::new(), Text::new());
    for i in 0..200 {
        let n = r.randint(0, 30) as usize;
        let s = r.chars("ACGTNacgu", n);
        let l = conv_label(r, i);
        let q = |r: &mut Rng, off: i64, lo: i64, hi: i64| -> String {
            (0..n)
                .map(|_| (off + r.randint(lo, hi)) as u8 as char)
                .collect()
        };
        let (a, b, c) = (q(r, 33, 0, 93), q(r, 64, 0, 62), q(r, 64, -5, 62));
        q33.push(&fastq(&l, &s, &a));
        q64.push(&fastq(&l, &s, &b));
        sol.push(&fastq(&l, &s, &c));
        if n > 0 {
            fa.push(&fasta(&l.replace('\t', " "), &s));
        }
    }
    w(t, "q33.fq", &q33.0)?;
    w(t, "q64.fq", &q64.0)?;
    w(t, "sol.fq", &sol.0)?;
    w(t, "a.fa", &fa.0)?;
    let bad = spoil(r, &q33.0, 4);
    w(t, "bad.fq", &bad)?;
    let bad = spoil(r, &fa.0, 2);
    w(t, "bad.fa", &bad)
}

fn pad8(mut b: Vec<u8>) -> Vec<u8> {
    while !b.len().is_multiple_of(8) {
        b.push(0);
    }
    b
}

pub fn sff(t: &Path, r: &mut Rng) -> io::Result<()> {
    let flows: usize = r.choice(&[4, 100, 400]);
    let read = |r: &mut Rng| -> Vec<u8> {
        let n = r.randint(0, 60) as usize;
        let reps = r.randint(1, 3) as usize;
        let name = format!("IJK{}", r.randint(0, 1_000_000)).repeat(reps);
        let bases = r.chars("ACGTN", n);
        let quals: Vec<u8> = (0..n).map(|_| r.randint(0, 60) as u8).collect();
        let clip = |r: &mut Rng| -> u16 {
            let x = r.randint(0, n as i64) as u16;
            r.choice(&[0, 0, x])
        };
        let mut rh = Vec::new();
        rh.extend((((16 + name.len()).div_ceil(8) * 8) as u16).to_be_bytes());
        rh.extend((name.len() as u16).to_be_bytes());
        rh.extend((n as u32).to_be_bytes());
        for _ in 0..4 {
            rh.extend(clip(r).to_be_bytes());
        }
        rh.extend(name.bytes());
        let mut data = vec![1u8; 2 * flows];
        data.extend(vec![1u8; n]);
        data.extend(bases.bytes());
        data.extend(quals);
        let mut out = pad8(rh);
        out.extend(pad8(data));
        out
    };
    let nreads = r.randint(1, 40) as usize;
    let reads: Vec<Vec<u8>> = (0..nreads).map(|_| read(r)).collect();
    let head_len = (31 + flows + 4).div_ceil(8) * 8;
    let where_: &str = r.choice(&["none", "end", "middle"]);
    let mut body = Vec::new();
    let (mut off, mut ilen) = (0u64, 0u32);
    let k = if where_ == "middle" {
        r.randint(0, reads.len() as i64) as usize
    } else {
        reads.len()
    };
    let index = |r: &mut Rng, body: &mut Vec<u8>, off: &mut u64, ilen: &mut u32| {
        *off = (head_len + body.len()) as u64;
        let mut data = b".mft1.00".to_vec();
        data.extend(vec![0u8; r.randint(0, 20) as usize]);
        *ilen = data.len() as u32;
        body.extend(pad8(data));
    };
    for (i, x) in reads.iter().enumerate() {
        if where_ != "none" && i == k {
            index(r, &mut body, &mut off, &mut ilen);
        }
        body.extend(x);
    }
    if where_ != "none" && k == reads.len() {
        index(r, &mut body, &mut off, &mut ilen);
    }
    let mut h = Vec::new();
    h.extend(0x2e73_6666u32.to_be_bytes());
    h.extend(1u32.to_be_bytes());
    h.extend(off.to_be_bytes());
    h.extend(ilen.to_be_bytes());
    h.extend((reads.len() as u32).to_be_bytes());
    h.extend((head_len as u16).to_be_bytes());
    h.extend(4u16.to_be_bytes());
    h.extend((flows as u16).to_be_bytes());
    h.push(1);
    h.extend(&b"TACG".repeat(flows)[..flows]);
    h.extend(b"TCAG");
    let mut file = pad8(h);
    file.extend(body);
    wb(t, "r.sff", &file)
}

fn ann(r: &mut Rng, bad: bool) -> String {
    let k: &str = r.choice(&["size=", "ee=", "length=", "sample=", "x="]);
    let v: &str = match k {
        "size=" => {
            let mut vals = vec!["1", "7", "12", "", "x", "+3", "1.5"];
            if bad {
                vals.extend(["0", "9223372036854775808", "99999999999999999999"]);
            }
            r.choice(&vals)
        }
        "ee=" => r.choice(&["0.25", "1", "1.2.3", ".", "", "2x"]),
        _ => r.choice(&["5", "", "ab"]),
    };
    format!("{k}{v}")
}

fn relab_label(r: &mut Rng, i: usize, bad: bool) -> String {
    let first = r.choice(&[format!("r{i}"), String::new(), format!("r{i} d")]);
    let mut parts = vec![first];
    for _ in 0..r.randint(0, 4) {
        parts.push(ann(r, bad));
    }
    let mut l = if r.chance(0.8) {
        parts.join(";")
    } else {
        parts.join(" ")
    };
    if r.chance(0.2) {
        l += ";";
    }
    if r.chance(0.1) {
        l = format!(";{l}");
    }
    if r.chance(0.2) {
        l += r.choice(&[" desc;size=4", " a b", "\tt;"]);
    }
    l
}

pub fn relab(t: &Path, r: &mut Rng) -> io::Result<()> {
    let recs = |r: &mut Rng, fq: bool, offset: i64, top: i64, empty: bool, bad: bool| {
        let mut out = Text::new();
        for i in 0..300 {
            let n = r.randint(if empty { 0 } else { 1 }, 40) as usize;
            let s = r.chars("ACGTNacgu", n);
            let l = relab_label(r, i, bad);
            if fq {
                let q: String = (0..n)
                    .map(|_| (offset + r.randint(0, top)) as u8 as char)
                    .collect();
                out.push(&fastq(&l, &s, &q));
            } else {
                out.push(&fasta(&l.replace('\t', " "), &s));
            }
        }
        out.0
    };
    let q33 = recs(r, true, 33, 93, false, false);
    let q64 = recs(r, true, 64, 62, false, false);
    let fa = recs(r, false, 0, 0, false, false);
    let e = recs(r, true, 33, 93, true, false);
    let bad = recs(r, false, 0, 0, false, true);
    w(t, "q33.fq", &q33)?;
    w(t, "q64.fq", &q64)?;
    w(t, "a.fa", &fa)?;
    w(t, "e.fq", &e)?;
    w(t, "bad.fa", &bad)
}

/// Reads whose scores wander, from `top` down.
#[allow(clippy::too_many_arguments)]
fn wandering(
    r: &mut Rng,
    n_reads: usize,
    offset: i64,
    top: i64,
    fq: bool,
    start: i64,
    step: (i64, i64),
    alphabet: &str,
    tag: &str,
    label: impl Fn(&mut Rng, usize) -> String,
) -> String {
    let mut out = Text::new();
    for i in 0..n_reads {
        let n = r.randint(0, if alphabet == "ACGTN" { 60 } else { 40 }) as usize;
        let mut s: Vec<u8> = (0..n).map(|_| r.ch(alphabet)).collect();
        if alphabet != "ACGTN" && r.chance(0.3) {
            // soft mask a stretch
            let a = r.randint(0, n.saturating_sub(1) as i64) as usize;
            let b = r.randint(a as i64, n as i64) as usize;
            s[a..b].make_ascii_lowercase();
        }
        let mut q = r.randint(start, top);
        let qs: String = s
            .iter()
            .map(|_| {
                q = (q + r.randint(step.0, step.1)).clamp(0, top);
                (offset + q) as u8 as char
            })
            .collect();
        let lab = label(r, i) + tag;
        let s = String::from_utf8(s).unwrap();
        out.push(&if fq {
            fastq(&lab, &s, &qs)
        } else {
            fasta(&lab, &s)
        });
    }
    out.0
}

pub fn sqtrim(t: &Path, r: &mut Rng) -> io::Result<()> {
    let label = |r: &mut Rng, i: usize| format!("r{i}") + if r.chance(0.2) { " d" } else { "" };
    let q33 = wandering(r, 400, 33, 41, true, 20, (-6, 4), "ACGTN", "", label);
    let hi = wandering(r, 400, 33, 93, true, 20, (-6, 4), "ACGTN", "", label);
    let q64 = wandering(r, 400, 64, 41, true, 20, (-6, 4), "ACGTN", "", label);
    let fa = wandering(r, 400, 0, 41, false, 20, (-6, 4), "ACGTN", "", label);
    w(t, "q33.fq", &q33)?;
    w(t, "hi.fq", &hi)?;
    w(t, "q64.fq", &q64)?;
    w(t, "a.fa", &fa)?;
    let bad = spoil(r, &q33, 4);
    w(t, "bad.fq", &bad)
}

pub fn sqfilt(t: &Path, r: &mut Rng) -> io::Result<()> {
    let al = "ACGTACGTACGTNn";
    let label = |r: &mut Rng, i: usize| {
        let mut l = format!("r{i}");
        if r.chance(0.5) {
            l += &size(r, 1, 9);
        }
        l
    };
    // the tag goes before the abundance and description, as /1 and /2 do
    let tagged = |tag: &'static str| {
        move |r: &mut Rng, i: usize| {
            let mut l = format!("r{i}{tag}");
            if r.chance(0.5) {
                l += &size(r, 1, 9);
            }
            if r.chance(0.1) {
                l += " d";
            }
            l
        }
    };
    let desc = |r: &mut Rng, i: usize| {
        let mut l = label(r, i);
        if r.chance(0.1) {
            l += " d";
        }
        l
    };
    let q33 = wandering(r, 400, 33, 41, true, 15, (-6, 5), al, "", desc);
    let q64 = wandering(r, 400, 64, 41, true, 15, (-6, 5), al, "", desc);
    let hi = wandering(r, 400, 33, 93, true, 15, (-6, 5), al, "", desc);
    let fa = wandering(r, 400, 0, 41, false, 15, (-6, 5), al, "", desc);
    let f = wandering(r, 400, 33, 41, true, 15, (-6, 5), al, "", tagged("/1"));
    let rv = wandering(r, 400, 33, 41, true, 15, (-6, 5), al, "", tagged("/2"));
    let hard_src = wandering(r, 400, 0, 41, false, 15, (-6, 5), al, "", desc);
    let hard: String = hard_src
        .split_inclusive('\n')
        .map(|l| {
            if l.starts_with('>') {
                l.to_string()
            } else {
                l.chars()
                    .map(|c| if c.is_ascii_lowercase() { 'N' } else { c })
                    .collect()
            }
        })
        .collect();
    w(t, "q33.fq", &q33)?;
    w(t, "q64.fq", &q64)?;
    w(t, "hi.fq", &hi)?;
    w(t, "a.fa", &fa)?;
    w(t, "f.fq", &f)?;
    w(t, "r.fq", &rv)?;
    w(t, "h.fa", &hard)
}

pub fn rerep(t: &Path, r: &mut Rng) -> io::Result<()> {
    let (mut fa, mut fq) = (Text::new(), Text::new());
    for i in 0..300 {
        let n = r.randint(0, 30) as usize;
        let s = r.chars("ACGTNacgu", n);
        let mut parts = vec![format!("r{i}")];
        for _ in 0..r.randint(0, 3) {
            let sz = format!("size={}", r.randint(1, 6));
            parts.push(r.choice(&[
                sz,
                "size=".into(),
                "size=x".into(),
                "x=1".into(),
                "ee=0.5".into(),
            ]));
        }
        let mut l = parts.join(";");
        if r.chance(0.2) {
            l += r.choice(&[" desc", ";", "\tt;size=2"]);
        }
        fa.push(&fasta(&l.replace('\t', " "), &s));
        let q: String = (0..n)
            .map(|_| (33 + r.randint(0, 41)) as u8 as char)
            .collect();
        fq.push(&fastq(&l, &s, &q));
    }
    w(t, "a.fa", &fa.0)?;
    w(t, "q.fq", &fq.0)
}

pub fn revcomp(t: &Path, r: &mut Rng) -> io::Result<()> {
    let codes = "ACGTUNRYKMBVDHSWacgtunrykmbvdhsw";
    let (mut fa, mut fq) = (Text::new(), Text::new());
    for i in 0..400 {
        let n = r.randint(0, 60) as usize;
        let s = r.chars(codes, n);
        let mut lab = format!("r{i}");
        if r.chance(0.3) {
            lab += &size(r, 1, 9);
        }
        if r.chance(0.2) {
            lab += " d";
        }
        fa.push(&fasta(&lab, &s));
        let q: String = (0..n)
            .map(|_| (33 + r.randint(0, 41)) as u8 as char)
            .collect();
        fq.push(&fastq(&lab, &s, &q));
    }
    w(t, "a.fa", &fa.0)?;
    w(t, "q.fq", &fq.0)
}

pub fn sortby(t: &Path, r: &mut Rng) -> io::Result<()> {
    let (mut fa, mut fq) = (Text::new(), Text::new());
    for _ in 0..500 {
        // many equal lengths
        let n: usize = r.choice(&[0, 1, 2, 3, 3, 4, 5, 8, 8, 8, 13, 40]);
        let s = r.chars("ACGTNacgu", n);
        let rn = format!("r{}", r.randint(0, 5));
        let mut lab: String = r
            .choice(&[
                "a",
                "A",
                "b",
                "B",
                "ab",
                "a b",
                "a;x=1",
                "zz",
                "Z",
                rn.as_str(),
            ])
            .to_string();
        if r.chance(0.5) {
            lab += &format!(";size={}", r.choice(&[1u64, 1, 2, 2, 3, 10, 4294967296]));
        }
        if r.chance(0.1) {
            lab += " desc";
        }
        fa.push(&fasta(&lab, &s));
        let q: String = (0..n)
            .map(|_| (33 + r.randint(0, 41)) as u8 as char)
            .collect();
        fq.push(&fastq(&lab, &s, &q));
    }
    w(t, "a.fa", &fa.0)?;
    w(t, "q.fq", &fq.0)
}

pub fn getseq(t: &Path, r: &mut Rng) -> io::Result<()> {
    let tokens = [
        "seq1", "SEQ1", "seq10", "s1", "S1", "s11", "x", "ok", "2", "a b", "",
    ];
    let header = |r: &mut Rng| {
        let parts: Vec<&str> = (0..r.randint(1, 3)).map(|_| r.choice(&tokens)).collect();
        let mut h = parts.join(r.choice(&[";", "-", ".", "|", " "]));
        if r.chance(0.4) {
            h += &format!(";sample={}", r.choice(&["s1", "s11", "S1", "s1;x"]));
        }
        if r.chance(0.3) {
            h += &format!(";tag={}", r.choice(&["seq1", "x", "ok"]));
        }
        if r.chance(0.2) {
            h += r.choice(&[" desc", ";", ";;"]);
        }
        h
    };
    let (mut fa, mut fq) = (Vec::new(), Text::new());
    for _ in 0..400 {
        let h = header(r);
        let n = r.randint(1, 10) as usize;
        let s = r.chars("ACGT", n);
        fa.push(fasta(&h, &s));
        fq.push(&fastq(&h, &s, &"I".repeat(n)));
    }
    w(t, "a.fa", &fa.concat())?;
    w(t, "q.fq", &fq.0)?;
    let labels: Vec<&str> = (0..5)
        .map(|_| {
            r.choice(&[
                "seq1",
                "SEQ1;sample=s1",
                "x",
                "ok",
                "",
                "s1-2",
                "seq10;tag=x",
            ])
        })
        .collect();
    w(t, "labels.txt", &(labels.join("\n") + "\n"))?;
    let words: Vec<&str> = (0..5)
        .map(|_| r.choice(&["seq1", "s1", "ok", "2", "x", "a b", "s1;x", ""]))
        .collect();
    w(t, "words.txt", &(words.join("\n") + "\n"))?;
    let mut bad = fa.clone();
    for i in r.sample(bad.len(), 3) {
        let spoil = r.choice(&[";size=0\n", ";size=9223372036854775808\n"]);
        bad[i] = bad[i].replacen('\n', spoil, 1);
    }
    w(t, "bad.fa", &bad.concat())
}

pub fn subseq(t: &Path, r: &mut Rng) -> io::Result<()> {
    let reads = |r: &mut Rng, bad: bool| {
        let (mut fa, mut fq) = (Text::new(), Text::new());
        for i in 0..300 {
            let n = r.randint(0, 40) as usize;
            let s = r.chars("ACGTNacgu", n);
            let mut lab = format!("r{i}");
            if r.chance(0.3) {
                lab += &size(r, 1, 9);
            }
            if r.chance(0.2) {
                lab += " d";
            }
            if bad && r.chance(0.01) {
                lab += ";size=0";
            }
            fa.push(&fasta(&lab, &s));
            let q: String = (0..n)
                .map(|_| (33 + r.randint(0, 41)) as u8 as char)
                .collect();
            fq.push(&fastq(&lab, &s, &q));
        }
        (fa.0, fq.0)
    };
    let (a, q) = reads(r, false);
    w(t, "a.fa", &a)?;
    w(t, "q.fq", &q)?;
    let (a, q) = reads(r, true);
    w(t, "bad.fa", &a)?;
    w(t, "bad.fq", &q)
}

pub fn fqchar(t: &Path, r: &mut Rng) -> io::Result<()> {
    // total residues of 400, 800 or 1600 give percentages that end in 25 or 75
    let total = {
        let x = r.randint(1, 3000);
        r.choice(&[400, 800, 1600, x])
    };
    let (lo, hi) = r.choice(&[
        (33, 73),
        (33, 74),
        (59, 104),
        (64, 104),
        (66, 104),
        (35, 40),
    ]);
    let mut out = Text::new();
    let (mut left, mut i) = (total, 0);
    while left > 0 || r.chance(0.1) {
        let n = left.min(r.randint(0, 60));
        left -= n;
        let s = r.chars("AAAACGTNacgtnuRY", n as usize);
        let mut q: Vec<u8> = (0..n).map(|_| r.randint(lo, hi) as u8).collect();
        if !q.is_empty() && r.chance(0.4) {
            // a tail of one symbol
            let k = r.randint(1, q.len().min(8) as i64) as usize;
            let c = r.randint(lo, hi) as u8;
            let len = q.len();
            q[len - k..].fill(c);
        }
        out.push(&fastq(&format!("r{i}"), &s, &String::from_utf8(q).unwrap()));
        i += 1;
        if i > 2000 {
            break;
        }
    }
    w(t, "q.fq", &out.0)
}

pub fn fqstat(t: &Path, r: &mut Rng) -> io::Result<()> {
    // reads of varied lengths, some long, with scores drawn from a narrow or
    // a wide band so that every cutoff is crossed somewhere; q.fq is
    // Phred+33, s.fq Phred+64 down to -5
    let reads = |r: &mut Rng, offset: i64, low: i64, high: i64, name: &str| {
        let n: usize = r.choice(&[1, 2, 5, 50, 500]);
        let top: i64 = r.choice(&[5, 30, 150, 700]);
        let mut out = Text::new();
        for i in 0..n {
            let x = r.randint(0, top);
            let len = r.choice(&[0, 1, 2, top, x]) as usize;
            let lo = r.randint(low, high);
            let hi = r.randint(lo, high);
            let s = r.chars("ACGTN", len);
            let q: String = (0..len)
                .map(|_| (offset + r.randint(lo, hi)) as u8 as char)
                .collect();
            out.push(&fastq(&format!("r{i}"), &s, &q));
        }
        w(t, name, &out.0)
    };
    reads(r, 33, 0, 45, "q.fq")?;
    reads(r, 64, -5, 45, "s.fq")
}

pub fn shufl(t: &Path, r: &mut Rng) -> io::Result<()> {
    let n = {
        let x = r.randint(1, 5000) as usize;
        r.choice(&[1, 2, 3, 10, 100, 1000, x])
    };
    let (mut fa, mut fq) = (Text::new(), Text::new());
    for i in 0..n {
        let len = r.randint(0, 20) as usize;
        let s = r.chars("ACGT", len);
        let mut lab = format!("r{i}");
        if r.chance(0.2) {
            lab += &size(r, 1, 9);
        }
        fa.push(&fasta(&lab, &s));
        fq.push(&fastq(&lab, &s, &"I".repeat(len)));
    }
    w(t, "a.fa", &fa.0)?;
    w(t, "q.fq", &fq.0)?;
    let seeds = [
        r.randint(1, 1 << 31).to_string(),
        r.randint(i64::MIN, -1).to_string(),
        r.randint(1, i64::MAX).to_string(),
        "9223372036854775807".into(),
        "-9223372036854775808".into(),
        "1".into(),
    ];
    w(t, "seeds", &seeds.join(" "))
}

pub fn subsmp(t: &Path, r: &mut Rng) -> io::Result<()> {
    let n = {
        let x = r.randint(1, 3000);
        r.choice(&[1, 5, 50, 500, x])
    };
    let (mut fa, mut fq) = (Text::new(), Text::new());
    for i in 0..n {
        let len = r.randint(0, 15) as usize;
        let s = r.chars("ACGT", len);
        let mut lab = format!("r{i}");
        if r.chance(0.7) {
            let x = r.randint(1, 50);
            lab += &format!(";size={}", r.choice(&[1, 2, 3, 10, 100, x]));
        }
        if r.chance(0.1) {
            lab += " d";
        }
        fa.push(&fasta(&lab, &s));
        fq.push(&fastq(&lab, &s, &"I".repeat(len)));
    }
    w(t, "a.fa", &fa.0)?;
    w(t, "q.fq", &fq.0)?;
    let mut draws = Vec::new();
    for _ in 0..6 {
        let seed = {
            let a = r.randint(1, 1 << 31);
            let b = r.randint(i64::MIN, -1);
            r.choice(&[a, b])
        };
        let count = format!("--count {}", r.randint(1, n));
        let digits = r.randint(0, 3) as i32;
        let pct = r.uniform(0.1, 100.0);
        let scale = 10f64.powi(digits);
        let pct = format!("--pct {}", py_float((pct * scale).round() / scale));
        let size = r.choice(&[count, pct, format!("--count {n}"), "--pct 100".into()]);
        draws.push(format!("{seed}|{size}"));
    }
    w(t, "draws", &(draws.join("\n") + "\n"))
}

pub fn scrmbl(t: &Path, r: &mut Rng) -> io::Result<()> {
    let (mut fa, mut fq) = (Text::new(), Text::new());
    for i in 0..200 {
        let n = {
            let x = r.randint(0, 400) as usize;
            r.choice(&[0, 1, 2, 5, 9, 10, 11, 30, 100, x])
        };
        let alphabet = r.choice(&["ACGT", "AC", "A", "ACGTN", "acgtACGT", "ACGTRYKMN"]);
        let s = r.chars(alphabet, n);
        let mut lab = format!("r{i}");
        if r.chance(0.2) {
            lab += &size(r, 1, 5);
        }
        fa.push(&fasta(&lab, &s));
        let q: String = (0..n)
            .map(|_| (33 + r.randint(0, 41)) as u8 as char)
            .collect();
        fq.push(&fastq(&lab, &s, &q));
    }
    w(t, "a.fa", &fa.0)?;
    w(t, "q.fq", &fq.0)?;
    let seeds = [
        r.randint(1, 1 << 31),
        r.randint(i64::MIN, -1),
        r.randint(1, i64::MAX),
    ];
    w(
        t,
        "seeds",
        &seeds
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(" "),
    )
}

pub fn lcmask(t: &Path, r: &mut Rng) -> io::Result<()> {
    // random stretches, repeats of short units and homopolymers, placed so
    // that regions fall anywhere in and across the 64-residue windows
    let piece = |r: &mut Rng| -> String {
        let kind = r.random();
        if kind < 0.4 {
            let alphabet = r.choice(&["ACGT", "ACGT", "ACGTN", "acgtACGT", "ACGTURYKMN", "AC"]);
            let n = r.randint(1, 80) as usize;
            return r.chars(alphabet, n);
        }
        if kind < 0.8 {
            let k = r.randint(1, 6) as usize;
            let unit = r.chars("ACGTacgtN", k);
            return unit.repeat(r.randint(1, 20) as usize);
        }
        let c = r.ch("ACGTNacgtn") as char;
        c.to_string().repeat(r.randint(1, 40) as usize)
    };
    let (mut fa, mut fq) = (Text::new(), Text::new());
    for i in 0..300 {
        let k: usize = r.choice(&[0, 1, 2, 3, 5, 10]);
        let mut s: String = (0..k).map(|_| piece(r)).collect();
        if r.chance(0.1) {
            let cut = r.randint(0, 12) as usize;
            s.truncate(cut);
        }
        let mut lab = format!("r{i}");
        if r.chance(0.2) {
            lab += &size(r, 1, 5);
        }
        if r.chance(0.1) {
            lab += " desc";
        }
        fa.push(&fasta(&lab, &s));
        let q: String = s
            .bytes()
            .map(|_| (33 + r.randint(0, 41)) as u8 as char)
            .collect();
        fq.push(&fastq(&lab, &s, &q));
    }
    w(t, "a.fa", &fa.0)?;
    w(t, "q.fq", &fq.0)
}

type Rec = (String, String);

pub fn syncpa(t: &Path, r: &mut Rng) -> io::Result<()> {
    // pairs named in several styles (/1 /2, _1 _2, .1 .2, Casava, none), the
    // reverse file shuffled, mates dropped on either side, and now and then
    // a duplicate label, a blank-only header or an empty file
    let styles = [
        ("/1", "/2"),
        ("_1", "_2"),
        (".1", ".2"),
        (" 1:N:0:1", " 2:N:0:1"),
        ("\t1", "\t2"),
        ("", ""),
        ("/1 x", "/2 y"),
        ("/3", "/3"),
    ];
    let n: usize = r.choice(&[0, 1, 5, 40, 300]);
    let (mut fwd, mut rev): (Vec<Rec>, Vec<Rec>) = (vec![], vec![]);
    for i in 0..n {
        let s = r.choice(&styles);
        let sized = format!("r{i};size={}", r.randint(0, 3));
        let base = if r.chance(0.02) {
            r.choice(&[format!("r{i}"), format!("r{i}"), sized, String::new()])
        } else {
            r.choice(&[format!("r{i}"), sized])
        };
        let seq = |r: &mut Rng| {
            let k = r.randint(0, 30) as usize;
            r.chars("ACGTacgtN", k)
        };
        if r.chance(0.85) {
            fwd.push((base.clone() + s.0, seq(r)));
        }
        if r.chance(0.85) {
            rev.push((base + s.1, seq(r)));
        }
    }
    if !fwd.is_empty() && r.chance(0.1) {
        let at = r.randrange(fwd.len() + 1);
        let x = r.choice(&fwd);
        fwd.insert(at, x);
    }
    if !rev.is_empty() && r.chance(0.1) {
        let at = r.randrange(rev.len() + 1);
        let x = r.choice(&rev);
        rev.insert(at, x);
    }
    r.shuffle(&mut rev);
    let fa = |rs: &[(String, String)]| rs.iter().map(|(h, s)| fasta(h, s)).collect::<String>();
    w(t, "f.fa", &fa(&fwd))?;
    w(t, "r.fa", &fa(&rev))?;
    let mut fq = |rs: &[(String, String)]| {
        rs.iter()
            .map(|(h, s)| {
                let q: String = s
                    .bytes()
                    .map(|_| (33 + r.randint(0, 40)) as u8 as char)
                    .collect();
                fastq(h, s, &q)
            })
            .collect::<String>()
    };
    let (a, b) = (fq(&fwd), fq(&rev));
    w(t, "f.fq", &a)?;
    w(t, "r.fq", &b)
}

pub fn joinpa(t: &Path, r: &mut Rng) -> io::Result<()> {
    // FASTQ pairs of any lengths in mixed-case IUPAC with U, headers with
    // descriptions, tabs and abundances (zero now and then, in either file),
    // and now and then one file longer than the other, or empty
    let n: usize = r.choice(&[0, 1, 5, 40, 300]);
    let header = |r: &mut Rng, i: usize| {
        let sized = format!("r{i};size={}", r.randint(1, 9));
        let h = r.choice(&[
            format!("r{i}"),
            format!("r{i}/1"),
            format!("r{i} x y"),
            format!("r{i}\tz"),
            sized,
            format!("r{i};size=x"),
        ]);
        if r.chance(0.01) { h + ";size=0" } else { h }
    };
    let read = |r: &mut Rng, i: usize| {
        let len = {
            let a = r.randint(0, 40) as usize;
            let b = r.randint(0, 300) as usize;
            r.choice(&[0, 1, a, b])
        };
        let s = r.chars("ACGTUNRYKMBVDHSWacgtunrykmbvdhsw", len);
        let h = header(r, i);
        let q: String = (0..len)
            .map(|_| (33 + r.randint(0, 60)) as u8 as char)
            .collect();
        fastq(&h, &s, &q)
    };
    let mut fwd: Vec<String> = (0..n).map(|i| read(r, i)).collect();
    let mut rev: Vec<String> = (0..n).map(|i| read(r, i)).collect();
    let x = r.random();
    if x < 0.05 && !fwd.is_empty() {
        fwd.pop();
    } else if x < 0.1 && !rev.is_empty() {
        rev.pop();
    } else if x < 0.15 {
        fwd.push(read(r, n));
    } else if x < 0.2 {
        rev.push(read(r, n));
    }
    w(t, "f.fq", &fwd.concat())?;
    w(t, "r.fq", &rev.concat())
}

pub fn mergpa(t: &Path, r: &mut Rng) -> io::Result<()> {
    // pairs read from both ends of random fragments: overlaps from none to
    // whole reads, staggered pairs, substitutions, Ns, lower case and other
    // IUPAC codes, tandem repeats, empty reads, scores from low to high, now
    // and then a zero abundance or one file longer than the other
    let frag = |r: &mut Rng, n: usize| -> String {
        if r.chance(0.1) {
            let k = r.randint(1, 6) as usize;
            let u = r.chars("ACGT", k);
            return u.repeat(n / u.len() + 1)[..n].to_string();
        }
        r.chars("ACGT", n)
    };
    let mutate = |r: &mut Rng, s: &str| -> String {
        let mut s: Vec<u8> = s.bytes().collect();
        for _ in 0..r.choice(&[0, 0, 1, 2, 5, 15]) {
            if !s.is_empty() {
                let p = r.randrange(s.len());
                s[p] = r.ch("ACGTNNRYacgtu");
            }
        }
        String::from_utf8(s).unwrap()
    };
    let qual = |r: &mut Rng, n: usize| -> String {
        let (a, b) = (r.randint(0, 45), r.randint(0, 45));
        let (lo, hi) = (a.min(b), a.max(b));
        (0..n)
            .map(|_| (33 + r.randint(lo, hi)) as u8 as char)
            .collect()
    };
    let n: usize = r.choice(&[0, 1, 5, 40, 200]);
    let (mut fwd, mut rev) = (Vec::new(), Vec::new());
    for i in 0..n {
        let l = {
            let a = r.randint(0, 30);
            let b = r.randint(20, 150);
            let c = r.randint(100, 300);
            r.choice(&[a, b, c])
        } as usize;
        let mut a = r.randint(0, l as i64) as usize;
        let mut b = r.randint(0, l as i64) as usize;
        if r.chance(0.6) {
            a = r.randint((l / 2) as i64, l as i64) as usize;
            b = r.randint((l / 2) as i64, l as i64) as usize;
        }
        let whole = frag(r, l + 20);
        let f = if r.chance(0.9) {
            whole[..a].to_string()
        } else {
            let k = r.randint(0, 10) as usize;
            whole[k..][..a.min(whole.len() - k)].to_string()
        };
        let g = if r.chance(0.9) {
            rc(&whole[l - b..l])
        } else {
            let back = r.randint(0, 15) as usize;
            rc(&whole[(l as i64 - b as i64 - back as i64).max(0) as usize..l])
        };
        let mut f = mutate(r, &f);
        let g = mutate(r, &g);
        if r.chance(0.2) {
            f = f.to_lowercase();
        }
        let mut h = format!("p{i}");
        if r.chance(0.2) {
            h += r.choice(&[";size=3", " x", "/1"]);
        }
        if r.chance(0.001) {
            h += ";size=0";
        }
        let (fq_, gq) = (qual(r, f.len()), qual(r, g.len()));
        fwd.push(fastq(&h, &f, &fq_));
        rev.push(fastq(&format!("p{i}/2"), &g, &gq));
    }
    let x = r.random();
    if x < 0.03 && !fwd.is_empty() {
        fwd.pop();
    } else if x < 0.06 && !rev.is_empty() {
        rev.pop();
    }
    w(t, "f.fq", &fwd.concat())?;
    w(t, "r.fq", &rev.concat())
}

pub fn redig(t: &Path, r: &mut Rng) -> io::Result<()> {
    // sequences in mixed case with IUPAC codes, now and then planted with a
    // restriction site or a short repeat, empty, or shorter than a pattern;
    // and a random pattern of 1 to 8 codes with ^ and _ anywhere
    let codes = "ACGTUNRYKMBVDHSW";
    let seq = |r: &mut Rng| -> String {
        let n = {
            let a = r.randint(0, 40) as usize;
            let b = r.randint(20, 300) as usize;
            r.choice(&[0, 1, 3, a, b])
        };
        let mut s: Vec<u8> = (0..n)
            .map(|_| {
                let al = if r.chance(0.9) { "ACGT" } else { codes };
                r.ch(al)
            })
            .collect();
        for _ in 0..r.randint(0, 4) {
            if !s.is_empty() {
                let p = r.randrange(s.len());
                let site: &str = r.choice(&["GAATTC", "GGATCC", "AAGCTT", "AAAA", "GCGC", "CCGG"]);
                s.splice(p..p, site.bytes());
            }
        }
        let s = String::from_utf8(s).unwrap();
        if r.chance(0.1) { s.to_lowercase() } else { s }
    };
    let n: usize = r.choice(&[0, 1, 5, 50, 300]);
    let mut hs: Vec<String> = (0..n)
        .map(|i| {
            let mut h = format!("s{i}");
            if r.chance(0.2) {
                h += r.choice(&[";size=2", " desc", "\tx"]);
            }
            h
        })
        .collect();
    if !hs.is_empty() && r.chance(0.02) {
        let k = r.randrange(hs.len());
        hs[k] += ";size=0";
    }
    let text: String = hs
        .iter()
        .map(|h| {
            let s = seq(r);
            format!(">{h}\n{s}\n")
        })
        .collect();
    w(t, "in.fa", &text)?;
    let known = [
        "G^AATT_C", "G^GATC_C", "A^AGCT_T", "CCG^G_", "^GAATTC_", "_GAATTC^", "G_GT^T", "A^A_A",
        "N^N_",
    ];
    let p = if r.chance(0.5) {
        r.choice(&known).to_string()
    } else {
        let all = format!("{codes}{}", codes.to_lowercase());
        let k = r.randint(1, 8) as usize;
        let mut b: Vec<char> = (0..k).map(|_| r.ch(&all) as char).collect();
        let at = r.randint(0, b.len() as i64) as usize;
        b.insert(at, '^');
        let at = r.randint(0, b.len() as i64) as usize;
        b.insert(at, '_');
        b.into_iter().collect()
    };
    w(t, "pattern", &p)
}
