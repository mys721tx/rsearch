//! Generators for the analysis tools' cases: databases and queries related
//! by substitutions, indels, repeats and chimeras, UDB files, taxonomies,
//! and the drawn options of each search.

use std::io::{self, Write};
use std::path::Path;
use std::process::{Command, Stdio};

use super::{
    edit, fasta_all, fastq, fastq_i, flags, lower_span, mutate, options, rc, rnd, slice, sprinkle,
    substitute, w, wb,
};
use crate::rng::Rng;
use crate::translate::py_float;

const MASKS: &[&str] = &[
    "dust", "dust", "none", "soft", "dusthard", "softhard", "nonehard",
];

/// The search criteria vsearch's searches share.
const CRITERIA: &[(&str, &[&str])] = &[
    ("maxdiffs", &["0", "2", "10"]),
    ("maxsubs", &["0", "3"]),
    ("maxgaps", &["0", "1", "3"]),
    ("mincols", &["10", "50"]),
    ("query_cov", &["0.5", "0.9"]),
    ("target_cov", &["0.5", "0.9"]),
    ("maxid", &["0.9", "0.99"]),
    ("mid", &["90", "99"]),
    ("idprefix", &["3", "8"]),
    ("idsuffix", &["2", "6"]),
    ("maxqt", &["0.9", "1.1"]),
    ("minqt", &["0.9", "1.1"]),
    ("minsl", &["0.8", "0.95"]),
    ("maxsl", &["0.9", "0.99"]),
    ("maxqsize", &["1", "3"]),
    ("mintsize", &["1", "3"]),
    ("minsizeratio", &["0.5", "2"]),
    ("maxsizeratio", &["0.5", "2"]),
];

const GAPOPEN: &[&str] = &[
    "10I/1E",
    "20I/2E",
    "5Q/8T",
    "*E/20I",
    "3I/1LQ/2RT",
    "0E/20I",
    "*I/2E",
];
const GAPEXT: &[&str] = &["1I/1E", "2I/1E", "1Q/2T", "2I/3E", "1QL/2QR/4TL/8TR/10I"];

fn scores(r: &mut Rng, o: &mut Vec<String>, p: f64) {
    if r.chance(p) {
        let (m, x) = (r.choice(&["1", "3", "5"]), r.choice(&["-1", "-2", "-6"]));
        o.extend(["match".into(), m.into(), "mismatch".into(), x.into()]);
    }
}

fn gaps(r: &mut Rng, o: &mut Vec<String>, p: f64, open: &[&str], ext: &[&str]) {
    if r.chance(p) {
        let (a, b) = (r.choice(open), r.choice(ext));
        o.extend(["gapopen".into(), a.into(), "gapext".into(), b.into()]);
    }
}

fn params(t: &Path, o: &[String]) -> io::Result<()> {
    w(t, "params", &(o.join(" ") + "\n"))
}

/// A random quality string for `s` (scores 0 to 40).
fn qual(r: &mut Rng, s: &str) -> String {
    s.bytes().map(|_| r.randint(33, 73) as u8 as char).collect()
}

/// Queries as FASTQ (with random scores) or FASTA, three in ten FASTQ;
/// returns the format.
fn queries(
    t: &Path,
    r: &mut Rng,
    recs: &[(String, String)],
    random_scores: bool,
) -> io::Result<&'static str> {
    if r.chance(0.3) {
        let text = if random_scores {
            recs.iter()
                .map(|(h, s)| {
                    let q = qual(r, s);
                    fastq(h, s, &q)
                })
                .collect()
        } else {
            fastq_i(recs)
        };
        w(t, "q.fq", &text)?;
        Ok("fq")
    } else {
        w(t, "q.fa", &fasta_all(recs))?;
        Ok("fa")
    }
}

/// A zero abundance on a random record, now and then.
fn zero_size(r: &mut Rng, hs: &mut [String]) {
    if !hs.is_empty() && r.chance(0.03) {
        let k = r.randrange(hs.len());
        hs[k] += ";size=0";
    }
}

pub fn mkudb(t: &Path, r: &mut Rng, word: Option<u32>) -> io::Result<()> {
    // sequences in mixed case with IUPAC codes, U, tandem repeats and
    // homopolymers DUST masks, repeats of one fragment (words shared by many
    // sequences), empty or short ones, now and then one over 50000; headers
    // with descriptions and abundances (zero now and then); FASTA or FASTQ;
    // a word length, length filters and one of vsearch's masking choices
    let codes = "ACGTUNRYKMBVDHSW";
    let n = r.randint(5, 40) as usize;
    let shared = rnd(r, n);
    let piece = |r: &mut Rng| -> String {
        let k = r.random();
        if k < 0.4 {
            let n = r.randint(0, 80) as usize;
            return (0..n)
                .map(|_| {
                    let al = if r.chance(0.95) { "ACGT" } else { codes };
                    r.ch(al) as char
                })
                .collect();
        }
        if k < 0.55 {
            let n = r.randint(1, 6) as usize;
            let u = rnd(r, n);
            return u.repeat(r.randint(2, 20) as usize);
        }
        if k < 0.65 {
            let c = r.ch("ACGTN") as char;
            return c.to_string().repeat(r.randint(2, 40) as usize);
        }
        if k < 0.85 {
            return shared.clone();
        }
        let n = r.randint(1, 30) as usize;
        r.chars("acgtun", n)
    };
    let seq = |r: &mut Rng| -> String {
        let k = r.random();
        if k < 0.05 {
            return String::new();
        }
        if k < 0.12 {
            let n = r.randint(1, 12) as usize;
            return rnd(r, n);
        }
        if k < 0.13 {
            let n = r.randint(50001, 60000) as usize;
            return rnd(r, n);
        }
        let mut s: String = (0..r.randint(1, 6)).map(|_| piece(r)).collect();
        if r.chance(0.05) {
            s = s.to_lowercase();
        }
        if r.chance(0.05) {
            s = s.replace('T', "U").replace('t', "u");
        }
        s
    };
    let n: usize = if word.is_none() {
        r.choice(&[0, 1, 2, 5, 9, 50, 300])
    } else {
        r.choice(&[1, 2, 9, 50])
    };
    let mut hs: Vec<String> = (0..n)
        .map(|i| {
            let mut h = format!("s{i}");
            if r.chance(0.2) {
                h += r.choice(&[";size=2", " desc x", "\tt", ";size=7;"]);
            }
            h
        })
        .collect();
    zero_size(r, &mut hs);
    let recs: Vec<(String, String)> = hs.into_iter().map(|h| (h, seq(r))).collect();
    let fmt = if r.chance(0.3) {
        let text: String = recs
            .iter()
            .map(|(h, s)| {
                let q = qual(r, s);
                fastq(h, s, &q)
            })
            .collect();
        w(t, "in.fq", &text)?;
        "fq"
    } else {
        w(t, "in.fa", &fasta_all(&recs))?;
        "fa"
    };
    let k = r.weighted(
        &[3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        &[8., 8., 8., 8., 8., 20., 8., 8., 6., 4., 1.],
    );
    // a word length given, for the long words (files of 4^k words)
    let k = word.unwrap_or(k);
    let x = r.randint(1, 60);
    let minlen = r.choice(&[0, 1, 32, x]);
    let x = r.randint(1, 200);
    let maxlen = if minlen > 0 {
        r.choice(&[50000, 50000, x])
    } else {
        1000000
    };
    let mask = r.choice(MASKS);
    w(
        t,
        "params",
        &format!("{k} {minlen} {maxlen} {mask} {fmt}\n"),
    )
}

pub fn udbinfo(t: &Path, r: &mut Rng) -> io::Result<()> {
    // a UDB header, valid or with one field wrong, in a file of random size
    // around 4 bytes a sequence; now and then a header cut short, or random
    // bytes
    let bits = |r: &mut Rng| r.getrandbits(32) as u32;
    let mut h: Vec<u32> = (0..50)
        .map(|_| if r.chance(0.3) { bits(r) } else { 0 })
        .collect();
    h[0] = 0x5544_4246;
    h[2] = 32;
    h[4] = r.randint(3, 15) as u32;
    h[17] = 0x746e;
    h[49] = 0x5544_4266;
    let x = bits(r);
    h[5] = r.choice(&[1, x]);
    let x = r.randint(0, 200) as u32;
    h[6] = r.choice(&[100, x]);
    let x = bits(r);
    h[11] = r.choice(&[0, x]);
    let (a, b, c) = (r.randint(1, 60) as u32, r.randint(1, 1000) as u32, bits(r));
    h[13] = r.choice(&[1, a, b, c]);
    if r.chance(0.3) {
        let i = r.choice(&[0, 2, 4, 13, 17, 49]);
        let x = bits(r);
        h[i] = r.choice(&[0, 1, 2, 16, 31, 64, 0x5544_4266, 0x5544_4246, x]);
    }
    let n = 4 * h[13] as i64;
    let (a, b) = (r.randint(0, 9), r.randint(0, 1000));
    let size = r.choice(&[200, n, n - 1, n + a, b]);
    let mut data: Vec<u8> = h.iter().flat_map(|v| v.to_le_bytes()).collect();
    for _ in 0..(size.min(5000) - 200).max(0) {
        data.push(r.getrandbits(8) as u8);
    }
    if r.chance(0.1) {
        data.truncate(r.randint(0, 199) as usize);
    }
    if r.chance(0.03) {
        data = (0..r.randint(0, 400))
            .map(|_| r.getrandbits(8) as u8)
            .collect();
    }
    wb(t, "in.udb", &data)
}

/// Runs `mkudb --word-len k` on FASTA text and returns the UDB file.
fn run_mkudb(mkudb: &Path, k: u32, fa: &str) -> io::Result<Vec<u8>> {
    let mut child = Command::new(mkudb)
        .args(["--word-len", &k.to_string()])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()?;
    let mut stdin = child.stdin.take().unwrap();
    let out = std::thread::scope(|s| {
        s.spawn(move || stdin.write_all(fa.as_bytes()));
        child.wait_with_output()
    })?;
    if !out.status.success() {
        return Err(io::Error::other("mkudb failed"));
    }
    Ok(out.stdout)
}

pub fn udbstat(t: &Path, r: &mut Rng, mkudb: &Path) -> io::Result<()> {
    // a UDB file mkudb writes for random sequences (shared fragments,
    // repeats, Ns, lower case, empty ones) at word length 3 to 10; then,
    // half the time, one part of it spoiled: a word of a header, a count,
    // the word list, an offset, a length, cut short anywhere, or bytes after
    // its end
    let n = r.randint(5, 30) as usize;
    let shared = rnd(r, n);
    let seq = |r: &mut Rng| -> String {
        let k = r.random();
        if k < 0.05 {
            return String::new();
        }
        if k < 0.3 {
            let n = r.randint(0, 20) as usize;
            return shared.clone() + &rnd(r, n);
        }
        if k < 0.4 {
            let c = r.ch("ACGT") as char;
            return c.to_string().repeat(r.randint(1, 30) as usize);
        }
        let n = r.randint(1, 300) as usize;
        (0..n)
            .map(|_| {
                let al = if r.chance(0.9) { "ACGTN" } else { "acgt" };
                r.ch(al) as char
            })
            .collect()
    };
    let seqs: usize = r.choice(&[1, 2, 3, 8, 20, 100, 400]);
    let fa: String = (0..seqs)
        .map(|i| {
            let d = r.choice(&["", " x", ";size=3"]);
            let s = seq(r);
            format!(">s{i}{d}\n{s}\n")
        })
        .collect();
    let k: u32 = r.choice(&[3, 3, 4, 5, 6, 7, 8, 8, 9, 10]);
    let mut d = run_mkudb(mkudb, k, &fa)?;
    let words = 4usize.pow(k);
    let get = |d: &[u8], at: usize| u32::from_le_bytes(d[at..at + 4].try_into().unwrap());
    let put = |d: &mut [u8], at: usize, v: u32| d[at..at + 4].copy_from_slice(&v.to_le_bytes());
    let total: usize = (0..words).map(|i| get(&d, 200 + 4 * i) as usize).sum();
    let magic = 200 + 4 * words;
    let list = 204 + 4 * words;
    let second = list + 4 * total;
    let offsets = second + 32;
    let hc = u64::from_le_bytes(d[second + 20..second + 28].try_into().unwrap());
    let headers = offsets + 4 * seqs;
    let lengths = headers + hc as usize;
    let sections = [
        200,
        204,
        magic,
        list,
        second,
        offsets,
        headers,
        lengths,
        lengths + 4 * seqs,
    ];
    let seqs32 = seqs as u32;
    let mut what = "none";
    if r.chance(0.5) {
        what = r.choice(&[
            "header", "count", "magic", "entry", "second", "offset", "length", "cut", "append",
        ]);
        match what {
            "header" => {
                let i = r.choice(&[0, 2, 4, 5, 6, 11, 13, 17, 49]);
                let (v, l, x) = (
                    get(&d, 4 * i),
                    (d.len() / 4) as u32,
                    r.getrandbits(32) as u32,
                );
                let v = r.choice(&[
                    0,
                    1,
                    2,
                    16,
                    v.wrapping_add(1),
                    v.wrapping_sub(1),
                    l,
                    l + 1,
                    x,
                ]);
                put(&mut d, 4 * i, v);
            }
            "count" => {
                let i = r.randrange(words);
                let v = get(&d, 200 + 4 * i);
                let v = r.choice(&[seqs32, seqs32 + 1, v + 1, v.saturating_sub(1), 0]);
                put(&mut d, 200 + 4 * i, v);
            }
            "magic" => {
                let v = r.getrandbits(32) as u32;
                put(&mut d, magic, v);
            }
            "entry" if total > 0 => {
                let at = list + 4 * r.randrange(total);
                let x = r.getrandbits(32) as u32;
                let v = r.choice(&[seqs32, seqs32 - 1, 0, x]);
                put(&mut d, at, v);
            }
            "second" => {
                let at = second + 4 * r.randrange(8);
                let v = get(&d, at).wrapping_add(r.choice(&[1u32, u32::MAX, 7]));
                put(&mut d, at, v);
            }
            "offset" => {
                let at = offsets + 4 * r.randrange(seqs);
                let v = get(&d, at);
                let v = r.choice(&[
                    v.wrapping_add(1),
                    v.wrapping_sub(1),
                    0,
                    hc as u32,
                    (hc as u32).wrapping_sub(1),
                ]);
                put(&mut d, at, v);
            }
            "length" => {
                let at = lengths + 4 * r.randrange(seqs);
                let v = get(&d, at).wrapping_add(r.choice(&[1u32, u32::MAX]));
                put(&mut d, at, v);
            }
            "cut" => {
                let (a, b) = (r.randrange(d.len()), r.choice(&sections));
                let at = r.choice(&[a, b, d.len() - 1]);
                d.truncate(at);
            }
            "append" => {
                for _ in 0..r.randint(1, 8) {
                    d.push(r.getrandbits(8) as u8);
                }
            }
            _ => {}
        }
    }
    wb(t, "in.udb", &d)?;
    w(t, "what", &format!("{what}\n"))
}

pub fn dbornt(t: &Path, r: &mut Rng) -> io::Result<()> {
    // a database of random sequences, tandem repeats and homopolymers in
    // mixed case with Ns, IUPAC codes and U, empty ones and now and then one
    // over 50000; queries cut from it on either strand, with substitutions,
    // lower case, chimeras of two strands, random or empty ones; abundances
    // (zero now and then) on both sides; FASTA or FASTQ queries; a word
    // length, a masking choice for each side, and a FASTA or UDB database
    let dbseq = |r: &mut Rng| -> String {
        let k = r.random();
        if k < 0.03 {
            return String::new();
        }
        if k < 0.04 {
            let n = r.randint(50001, 52000) as usize;
            return rnd(r, n);
        }
        let mut s = String::new();
        for _ in 0..r.randint(1, 5) {
            s += &match r.randrange(4) {
                0 => {
                    let n = r.randint(5, 200) as usize;
                    rnd(r, n)
                }
                1 => {
                    let n = r.randint(2, 5) as usize;
                    let u = rnd(r, n);
                    u.repeat(r.randint(3, 20) as usize)
                }
                2 => {
                    let c = r.ch("ACGT") as char;
                    c.to_string().repeat(r.randint(5, 30) as usize)
                }
                _ => {
                    let n = r.randint(1, 30) as usize;
                    r.chars("ACGTNRY", n)
                }
            };
        }
        if r.chance(0.1) {
            s = s.to_lowercase();
        } else if r.chance(0.2) {
            s = lower_span(r, &s);
        }
        if r.chance(0.05) {
            s = s.replace('T', "U");
        }
        s
    };
    let nd: usize = r.choice(&[0, 1, 3, 10, 50, 200]);
    let db: Vec<String> = (0..nd).map(|_| dbseq(r)).collect();
    let mut dh: Vec<String> = (0..nd)
        .map(|i| {
            let mut h = format!("d{i}");
            if r.chance(0.1) {
                h += r.choice(&[";size=4", " x"]);
            }
            h
        })
        .collect();
    zero_size(r, &mut dh);
    let good: Vec<&String> = db.iter().filter(|s| s.len() > 20).collect();
    let query = |r: &mut Rng| -> String {
        let k = r.random();
        if k < 0.1 || good.is_empty() {
            if r.chance(0.9) {
                let n = r.randint(0, 150) as usize;
                return rnd(r, n);
            }
            return String::new();
        }
        let s = r.choice(&good);
        let i = r.randrange(s.len());
        let mut q = slice(s, i, r.randint(5, 300) as usize).to_string();
        if k < 0.2 && good.len() > 1 {
            let t = r.choice(&good);
            let j = r.randrange(t.len());
            q += &rc(slice(t, j, r.randint(5, 200) as usize));
        }
        if r.chance(0.5) {
            let k = r.randint(0, (q.len() / 10 + 1) as i64) as usize;
            let mut b = q.into_bytes();
            for _ in 0..k {
                let p = r.randrange(b.len());
                b[p] = r.ch("ACGTN");
            }
            q = String::from_utf8(b).unwrap();
        }
        if r.chance(0.5) {
            q = rc(&q);
        }
        if r.chance(0.1) {
            q = q.to_lowercase();
        } else if r.chance(0.1) {
            let i = r.randrange(q.len() + 1);
            q = q[..i].to_lowercase() + &q[i..];
        }
        q
    };
    let nq: usize = r.choice(&[0, 1, 5, 30, 200]);
    let qs: Vec<String> = (0..nq).map(|_| query(r)).collect();
    let mut qh: Vec<String> = (0..nq)
        .map(|i| {
            let mut h = format!("q{i}");
            if r.chance(0.15) {
                h += r.choice(&[";size=2", " desc", "\tt"]);
            }
            h
        })
        .collect();
    zero_size(r, &mut qh);
    let dbrecs: Vec<(String, String)> = dh.into_iter().zip(db.iter().cloned()).collect();
    w(t, "db.fa", &fasta_all(&dbrecs))?;
    let qrecs: Vec<(String, String)> = qh.into_iter().zip(qs).collect();
    let fmt = queries(t, r, &qrecs, true)?;
    let k = r.weighted(
        &[3, 4, 5, 6, 8, 10, 11, 12, 13],
        &[5., 5., 5., 8., 10., 8., 8., 40., 3.],
    );
    let dbmask = r.choice(MASKS);
    let qmask = r.choice(&["dust", "dust", "soft", "none"]);
    let dbkind = r.choice(&["fasta", "fasta", "fasta", "udb"]);
    w(
        t,
        "params",
        &format!("{k} {dbmask} {qmask} {dbkind} {fmt}\n"),
    )
}

pub fn sintax(t: &Path, r: &mut Rng) -> io::Result<()> {
    // a taxonomy of templates in a random tree; a database of copies of
    // them, mutated, truncated (ties broken by length), duplicated, with
    // ;tax= in many spellings (missing levels, upper-case letters, empty
    // names, other annotations around it, none at all), stretches in lower
    // case, short ones and now and then one over 50000; queries cut from the
    // templates on either strand, mutated, chimeric, random, short, in lower
    // case; zero abundances now and then; and the options
    let mutate_at = |r: &mut Rng, s: &str, rate: f64| -> String {
        s.chars()
            .map(|c| {
                if r.chance(rate) {
                    r.ch("ACGT") as char
                } else {
                    c
                }
            })
            .collect()
    };
    let nt = r.randint(1, 8) as usize;
    let len: usize = r.choice(&[80, 200, 500]);
    let tem: Vec<String> = (0..nt).map(|_| rnd(r, len)).collect();
    // names by level, shared along random tree paths
    let levels = b"dkpcofgst";
    let mut paths: Vec<Vec<String>> = Vec::new();
    for i in 0..nt {
        let p = if i > 0 && r.chance(0.6) {
            let mut p = paths[r.randrange(i)].clone();
            let cut = r.randint(1, 8) as usize;
            p.truncate(cut);
            p.extend((cut..9).map(|j| format!("{}{i}_{j}", levels[j] as char)));
            p
        } else {
            (0..9)
                .map(|j| format!("{}{i}_{j}", levels[j] as char))
                .collect()
        };
        paths.push(p);
    }
    let taxstr = |r: &mut Rng, i: usize| -> String {
        let mut parts = Vec::new();
        for (j, name) in paths[i].iter().enumerate() {
            if r.chance(0.25) {
                continue;
            }
            let l = levels[j] as char;
            let letter = if r.chance(0.05) {
                l.to_ascii_uppercase()
            } else {
                l
            };
            if r.chance(0.03) {
                parts.push(format!("{letter}:"));
                continue;
            }
            parts.push(format!("{letter}:{name}"));
        }
        if r.chance(0.05) {
            let at = r.randrange(parts.len() + 1);
            parts.insert(at, "z:junk".into());
        }
        parts.join(",")
    };
    let mut db: Vec<(String, String)> = Vec::new();
    for i in 0..r.choice(&[0, 1, 5, 20, 60]) {
        let tx = r.randrange(nt);
        let rate = r.choice(&[0.0, 0.0, 0.01, 0.05, 0.15]);
        let mut s = mutate_at(r, &tem[tx], rate);
        let k = r.random();
        if k < 0.1 {
            s.truncate(r.randint(1, s.len() as i64) as usize);
        }
        if k > 0.97 {
            let n = r.randint(50001, 51000) as usize;
            s = rnd(r, n);
        }
        if r.chance(0.1) {
            s = lower_span(r, &s);
        }
        let mut h = format!("r{i}");
        let k = r.random();
        if k < 0.8 {
            h += &format!(";tax={};", taxstr(r, tx));
        } else if k < 0.9 {
            h += &format!(";x=1;tax={}", taxstr(r, tx));
        } else if k < 0.95 {
            h += &format!(" desc;tax={};size=3", taxstr(r, tx));
        }
        let dup = r.chance(0.1);
        db.push((h.clone(), s.clone()));
        if dup {
            db.push((h + "dup", s));
        }
    }
    if !db.is_empty() && r.chance(0.03) {
        let k = r.randrange(db.len());
        db[k] = (db[0].0.clone() + ";size=0", db[0].1.clone());
    }
    let mut qs = Vec::new();
    for i in 0..r.choice(&[0, 1, 10, 40]) {
        let k = r.random();
        let mut s;
        if k < 0.1 {
            let n = r.randint(0, 200) as usize;
            s = rnd(r, n);
        } else if k < 0.15 {
            let n = r.randint(0, 40) as usize;
            s = rnd(r, n);
        } else {
            let tm = &tem[r.randrange(nt)];
            let a = r.randrange(tm.len());
            let cut = slice(tm, a, r.randint(20, 300) as usize).to_string();
            let rate = r.choice(&[0.0, 0.02, 0.1]);
            s = mutate_at(r, &cut, rate);
            if k < 0.25 {
                let t2 = &tem[r.randrange(nt)];
                let h = s.len() / 2;
                s = s[..h].to_string() + &t2[..h];
            }
            if r.chance(0.4) {
                s = rc(&s);
            }
            if r.chance(0.1) {
                s = lower_span(r, &s);
            }
        }
        let mut h = format!("q{i}");
        if r.chance(0.1) {
            h += r.choice(&[" d", ";size=2", "\tt"]);
        }
        qs.push((h, s));
    }
    if !qs.is_empty() && r.chance(0.03) {
        let k = r.randrange(qs.len());
        qs[k] = ("z;size=0".into(), "ACGT".into());
    }
    let fmt = queries(t, r, &qs, false)?;
    w(t, "db.fa", &fasta_all(&db))?;
    let k = r.weighted(
        &[3, 4, 5, 6, 7, 8, 9, 10, 12, 13],
        &[2., 3., 4., 6., 6., 30., 6., 6., 6., 2.],
    );
    let strand = r.choice(&["plus", "plus", "both"]);
    let x = py_float((r.random() * 100.0).round() / 100.0);
    let cutoff = r.choice(&[
        "0".to_string(),
        "0".into(),
        "0.5".into(),
        "0.8".into(),
        "1".into(),
        x,
    ]);
    let rnd_ = r.choice(&[0, 0, 1]);
    let x = r.randint(1, 1_000_000_000);
    let seed = r.choice(&[1, 2, 3, -5, 1 << 40, x]);
    let dbmask = r.choice(&["dust", "soft", "none"]);
    let qmode = r.choice(&["mask", "mask", "clear"]);
    let dbkind = r.choice(&["fasta", "fasta", "udb"]);
    w(
        t,
        "params",
        &format!("{k} {strand} {cutoff} {rnd_} {seed} {dbmask} {qmode} {dbkind} {fmt}\n"),
    )
}

pub fn exsrch(t: &Path, r: &mut Rng) -> io::Result<()> {
    // a database with many copies of few sequences, in either case, with U,
    // IUPAC codes, reverse complements, empty ones and now and then one over
    // 50000; labels shared now and then, abundances (zero now and then);
    // queries copied from it, recased, reverse-complemented, mutated, random
    // or empty; and the options
    let base: Vec<String> = (0..r.randint(1, 8))
        .map(|_| {
            let n = r.randint(0, 30) as usize;
            let al = if r.chance(0.8) { "ACGT" } else { "ACGTNRY" };
            r.chars(al, n)
        })
        .collect();
    let variant = |r: &mut Rng, s: &str| -> String {
        let k = r.random();
        let mut s = if k < 0.2 {
            s.to_lowercase()
        } else if k < 0.3 {
            s.chars()
                .map(|c| {
                    if r.chance(0.5) {
                        c.to_ascii_lowercase()
                    } else {
                        c
                    }
                })
                .collect()
        } else {
            s.to_string()
        };
        if r.chance(0.1) {
            s = s.replace('T', "U").replace('t', "u");
        }
        if r.chance(0.15) {
            s = rc(&s);
        }
        s
    };
    let label = |r: &mut Rng, p: &str, i: usize| -> String {
        let mut h = if r.chance(0.8) {
            format!("{p}{i}")
        } else {
            format!("s{}", r.randrange(3))
        };
        if r.chance(0.3) {
            h += &format!(";size={}", r.choice(&[1, 2, 5, 10, 100]));
        }
        h
    };
    let mut db = Vec::new();
    for i in 0..r.choice(&[0, 1, 5, 20, 100]) {
        let k = r.random();
        let s = if k < 0.03 {
            String::new()
        } else if k < 0.04 {
            let n = r.randint(50001, 50100) as usize;
            rnd(r, n)
        } else if k < 0.8 {
            let b = r.choice(&base);
            variant(r, &b)
        } else {
            let n = r.randint(1, 30) as usize;
            rnd(r, n)
        };
        let h = label(r, "d", i);
        db.push((h, s));
    }
    if !db.is_empty() && r.chance(0.03) {
        let k = r.randrange(db.len());
        db[k] = ("z;size=0".into(), "ACGT".into());
    }
    let mut qs = Vec::new();
    for i in 0..r.choice(&[0, 1, 10, 50]) {
        let k = r.random();
        let mut s = if k < 0.6 {
            let b = r.choice(&base);
            variant(r, &b)
        } else if !db.is_empty() && k < 0.8 {
            r.choice(&db).1
        } else {
            let n = r.randint(0, 20) as usize;
            rnd(r, n)
        };
        if r.chance(0.1) && !s.is_empty() {
            let j = r.randrange(s.len());
            let c = r.ch("ACGT") as char;
            s = format!("{}{c}{}", &s[..j], &s[j + 1..]);
        }
        let h = label(r, "q", i);
        qs.push((h, s));
    }
    if !qs.is_empty() && r.chance(0.03) {
        let k = r.randrange(qs.len());
        qs[k] = ("z;size=0".into(), "ACGT".into());
    }
    let fmt = queries(t, r, &qs, false)?;
    w(t, "db.fa", &fasta_all(&db))?;
    let mut o = vec![fmt.to_string()];
    flags(r, &mut o, 0.4, &["strand"]);
    if r.chance(0.3) {
        o.extend(["maxhits".into(), r.randint(1, 4).to_string()]);
    }
    flags(r, &mut o, 0.3, &["self"]);
    let sizes: &[&str] = &["0", "1", "2", "5", "50"];
    let ratios: &[&str] = &["0.1", "0.5", "1", "2", "0.3333"];
    options(
        r,
        &mut o,
        0.15,
        &[
            ("maxqsize", sizes),
            ("mintsize", sizes),
            ("minsizeratio", ratios),
            ("maxsizeratio", ratios),
        ],
    );
    params(t, &o)
}

pub fn hitfmt(t: &Path, r: &mut Rng) -> io::Result<()> {
    // a database of related sequences (templates, copies with
    // substitutions, indels, IUPAC codes, lower case), and queries cut from
    // them on either strand with substitutions, insertions, deletions and
    // overhangs, or random; and the options of a vsearch search
    let len: usize = r.choice(&[60, 150, 300]);
    let tem: Vec<String> = (0..r.randint(1, 5)).map(|_| rnd(r, len)).collect();
    let mut db = Vec::new();
    for i in 0..r.randint(1, 30) {
        let tm = r.choice(&tem);
        let rate = r.choice(&[0.0, 0.02, 0.08]);
        let mut s = edit(r, &tm, rate, true);
        if r.chance(0.1) && !s.is_empty() {
            let a = r.randrange(s.len());
            let b = (a + 20).min(s.len());
            s = format!("{}{}{}", &s[..a], s[a..b].to_lowercase(), &s[b..]);
        }
        let d = if r.chance(0.1) { " x" } else { "" };
        db.push((format!("d{i}{d}"), s));
    }
    let mut qs = Vec::new();
    for i in 0..r.randint(1, 30) {
        let s = if r.chance(0.1) {
            let n = r.randint(32, 200) as usize;
            rnd(r, n)
        } else {
            let tm = r.choice(&tem);
            let a = r.randrange(tm.len() / 2);
            let cut = slice(&tm, a, r.randint(40, len as i64) as usize).to_string();
            let rate = r.choice(&[0.0, 0.03, 0.1]);
            let mut s = edit(r, &cut, rate, true);
            if r.chance(0.3) {
                let n = r.randint(0, 15) as usize;
                s = rnd(r, n) + &s;
            }
            if r.chance(0.3) {
                let n = r.randint(0, 15) as usize;
                s += &rnd(r, n);
            }
            if r.chance(0.4) {
                s = rc(&s);
            }
            s
        };
        let d = if r.chance(0.1) { ";size=2" } else { "" };
        qs.push((format!("q{i}{d}"), s));
    }
    w(t, "db.fa", &fasta_all(&db))?;
    w(t, "q.fa", &fasta_all(&qs))?;
    let mut o: Vec<String> = vec![
        "--id".into(),
        r.choice(&["0.5", "0.7", "0.8", "0.9"]).into(),
        "--maxaccepts".into(),
        r.choice(&["1", "3", "10"]).into(),
        "--strand".into(),
        r.choice(&["plus", "both"]).into(),
    ];
    if r.chance(0.3) {
        o.extend(["--iddef".into(), r.randint(0, 4).to_string()]);
    }
    flags(r, &mut o, 0.2, &["--n_mismatch"]);
    if r.chance(0.2) {
        let (m, x) = (r.choice(&["1", "2", "3"]), r.choice(&["-1", "-2", "-4"]));
        o.extend(["--match".into(), m.into(), "--mismatch".into(), x.into()]);
    }
    // opening penalties at least their extensions. vsearch scores the
    // right-end gaps of its SIMD alignments otherwise in ways a CIGAR alone
    // does not give
    if r.chance(0.2) {
        let (a, b) = r.choice(&[
            ("10I/1E", "1I/1E"),
            ("20I/2E", "2I/1E"),
            ("5Q/8T", "1Q/2T"),
            ("*E/20I", "2I/1E"),
            ("3I/2LQ/2RT", "1Q/2T"),
        ]);
        o.extend(["--gapopen".into(), a.into(), "--gapext".into(), b.into()]);
    }
    flags(r, &mut o, 0.3, &["--output_no_hits"]);
    flags(r, &mut o, 0.2, &["--top_hits_only"]);
    if r.chance(0.2) {
        o.extend(["--maxhits".into(), r.randint(1, 3).to_string()]);
    }
    let mask = r.choice(&["none", "none", "dust"]);
    let width = r.choice(&[64, 64, 30, 0, 7]);
    w(t, "params", &format!("{mask} {width} {}\n", o.join(" ")))
}

pub fn otutab(t: &Path, r: &mut Rng) -> io::Result<()> {
    // a database of related sequences whose labels name OTUs (otu=, or the
    // label, some shared) with taxonomies (tax=, empty, or none; the same
    // for an OTU's sequences); queries cut from them, mutated, on either
    // strand, labelled with samples (sample=, barcodelabel=, or a prefix)
    // and abundances; a search, exact or global with weak hits
    let tem: Vec<String> = (0..r.randint(1, 8))
        .map(|_| {
            let n = r.choice(&[40, 100]);
            rnd(r, n)
        })
        .collect();
    let otus: Vec<String> = (0..r.randint(1, 6))
        .map(|i| format!("O{}", (b'A' + i as u8) as char))
        .collect();
    let tax: Vec<Option<String>> = otus
        .iter()
        .map(|o| {
            r.choice(&[
                Some(String::new()),
                None,
                Some(format!("d:B,p:P{o}")),
                Some(format!("d:A,g:{o}")),
            ])
        })
        .collect();
    let mut db = Vec::new();
    for i in 0..r.randint(1, 20) {
        let oi = r.randrange(otus.len());
        let o = &otus[oi];
        let k = r.random();
        let mut h = if k < 0.4 {
            o.clone()
        } else if k < 0.7 {
            format!("s{i};otu={o}")
        } else {
            format!("{o};size={}", r.randint(1, 9))
        };
        if let Some(x) = &tax[oi] {
            h += &format!(";tax={x}");
        }
        let tm = r.choice(&tem);
        let n = r.choice(&[0, 0, 2, 10]);
        db.push((h, substitute(r, &tm, n)));
    }
    let smp: Vec<&str> = (0..4)
        .map(|_| r.choice(&["S1", "S2", "x_3", "Q9"]))
        .collect();
    let mut qs = Vec::new();
    for i in 0..r.randint(0, 40) {
        let tm = r.choice(&tem);
        let mut s = if r.chance(0.9) {
            let n = r.choice(&[0, 0, 1, 5, 15]);
            substitute(r, &tm, n)
        } else {
            rnd(r, tm.len())
        };
        if r.chance(0.3) {
            s = rc(&s);
        }
        let sm = r.choice(&smp);
        let k = r.random();
        let mut h = if k < 0.3 {
            format!("{sm}.{i}")
        } else if k < 0.5 {
            format!("r{i};sample={sm}")
        } else if k < 0.6 {
            format!("r{i};barcodelabel={sm};sample=zz")
        } else {
            format!("{sm}_{i}")
        };
        if r.chance(0.4) {
            h += &format!(";size={}", r.randint(1, 20));
        }
        qs.push((h, s));
    }
    w(t, "db.fa", &fasta_all(&db))?;
    w(t, "q.fa", &fasta_all(&qs))?;
    let mode = r.choice(&["exact", "global", "global"]);
    let mut o: Vec<String> = vec!["--strand".into(), r.choice(&["plus", "both"]).into()];
    if mode == "global" {
        let (id, acc) = (
            r.choice(&["0.8", "0.9", "0.97"]),
            r.choice(&["1", "2", "5"]),
        );
        o.extend(["--id".into(), id.into(), "--maxaccepts".into(), acc.into()]);
        if r.chance(0.5) {
            o.extend(["--weak_id".into(), r.choice(&["0.5", "0.7"]).into()]);
        }
    }
    w(t, "params", &format!("{mode} {}\n", o.join(" ")))
}

pub fn lcatax(t: &Path, r: &mut Rng) -> io::Result<()> {
    // a database of mutated copies of templates whose taxonomies share
    // levels on random tree paths (level letters in either case, empty
    // names, missing levels, no tax= at all); queries cut from the
    // templates; a global search with many hits, --top_hits_only now and
    // then, and a cutoff
    let levels = b"dkpcofgst";
    let mut paths: Vec<Vec<String>> = Vec::new();
    for i in 0..r.randint(1, 8) {
        let p = if !paths.is_empty() && r.chance(0.7) {
            let mut p = r.choice(&paths);
            let c = r.randint(0, 8) as usize;
            p.truncate(c);
            p.extend((c..9).map(|j| format!("{}{i}", levels[j] as char)));
            p
        } else {
            (0..9)
                .map(|j| format!("{}{i}", levels[j] as char))
                .collect()
        };
        paths.push(p);
    }
    let tax = |r: &mut Rng, p: &[String]| -> String {
        let mut parts = Vec::new();
        for (j, n) in p.iter().enumerate() {
            if r.chance(0.2) {
                continue;
            }
            let l = levels[j] as char;
            let l = if r.chance(0.05) {
                l.to_ascii_uppercase()
            } else {
                l
            };
            parts.push(if r.chance(0.03) {
                format!("{l}:")
            } else {
                format!("{l}:{n}")
            });
        }
        parts.join(",")
    };
    let tem: Vec<String> = (0..r.randint(1, 3)).map(|_| rnd(r, 100)).collect();
    let mut db = Vec::new();
    for i in 0..r.randint(1, 25) {
        let tm = r.randrange(tem.len());
        let mut h = format!("r{i}");
        if r.chance(0.9) {
            let p = r.choice(&paths);
            h += &format!(";tax={}", tax(r, &p));
        }
        let n = r.randint(0, 15) as usize;
        db.push((h, substitute(r, &tem[tm], n)));
    }
    let q: Vec<(String, String)> = (0..r.randint(1, 30))
        .map(|i| {
            let s = if r.chance(0.9) {
                let tm = r.choice(&tem);
                let n = r.randint(0, 8) as usize;
                substitute(r, &tm[5..95], n)
            } else {
                rnd(r, 90)
            };
            (format!("q{i}"), s)
        })
        .collect();
    w(t, "db.fa", &fasta_all(&db))?;
    w(t, "q.fa", &fasta_all(&q))?;
    let (id, acc) = (
        r.choice(&["0.7", "0.8", "0.9"]),
        r.choice(&["1", "3", "8", "20"]),
    );
    let mut o: Vec<String> = vec![
        "--id".into(),
        id.into(),
        "--maxaccepts".into(),
        acc.into(),
        "--maxrejects".into(),
        "64".into(),
    ];
    flags(r, &mut o, 0.3, &["--top_hits_only"]);
    if r.chance(0.3) {
        o.extend(["--weak_id".into(), "0.5".into()]);
    }
    let cutoff = r.choice(&["1.0", "1.0", "0.51", "0.6", "0.75", "0.9"]);
    w(t, "params", &format!("{cutoff} {}\n", o.join(" ")))
}

/// Templates of length `len`: random, or (one in five) a repeat of three
/// bases.
fn templates(r: &mut Rng, len: usize, n: usize) -> Vec<String> {
    (0..n)
        .map(|_| {
            if r.chance(0.8) {
                rnd(r, len)
            } else {
                rnd(r, 3).repeat(len / 3)
            }
        })
        .collect()
}

/// A piece of a template: up to a `part`th cut off each end, edited, and
/// now and then with random overhangs (each with probability `over`).
fn piece(
    r: &mut Rng,
    tem: &[String],
    part: usize,
    rates: &[f64],
    over: f64,
    over_len: i64,
) -> String {
    let tm = r.choice(tem);
    let a = r.randrange((tm.len() / part).max(1));
    let b = tm.len() - r.randrange((tm.len() / part).max(1));
    let cut = if b > a { &tm[a..b] } else { "" };
    let rate = r.choice(rates);
    let mut s = edit(r, cut, rate, false);
    if r.chance(over) {
        let n = r.randint(0, over_len) as usize;
        s = rnd(r, n) + &s;
    }
    if r.chance(over) {
        let n = r.randint(0, over_len) as usize;
        s += &rnd(r, n);
    }
    s
}

/// A label `{p}{i}`, or now and then (`1 - own`) one of two shared ones,
/// with an abundance from 1 to `max` with probability `sized`.
fn shared_label(r: &mut Rng, p: &str, i: usize, own: f64, sized: f64, max: i64) -> String {
    let mut h = if r.chance(own) {
        format!("{p}{i}")
    } else {
        format!("dup{}", r.randrange(2))
    };
    if r.chance(sized) {
        h += &format!(";size={}", r.randint(1, max));
    }
    h
}

/// A pair of long related sequences for the fallback aligner.
fn big(r: &mut Rng) -> (String, String) {
    let n = r.randint(5000, 5300) as usize;
    let b = rnd(r, n);
    let e = edit(r, &b, 0.01, false);
    (b, e)
}

pub fn allpr(t: &Path, r: &mut Rng) -> io::Result<()> {
    // sequences related by substitutions, insertions, deletions, overhangs
    // and repeats, in mixed case with IUPAC codes, empty or random ones, now
    // and then long ones for the fallback aligner; abundances and shared
    // labels; and options: identity or accept all, criteria, scores, gap
    // penalties (negative effective openings, infinite ones), N as mismatch
    let len: usize = r.choice(&[20, 60, 150, 400]);
    let n = r.randint(1, 3) as usize;
    let tem = templates(r, len, n);
    let mut seqs = Vec::new();
    for i in 0..r.randint(1, 20) as usize {
        let k = r.random();
        let mut s = if k < 0.05 {
            String::new()
        } else if k < 0.1 {
            let n = r.randint(1, len as i64) as usize;
            rnd(r, n)
        } else {
            piece(r, &tem, 4, &[0.0, 0.01, 0.05, 0.15, 0.3], 0.2, 20)
        };
        if r.chance(0.1) {
            s = sprinkle(r, &s, 0.05);
        }
        if r.chance(0.1) {
            s = s.to_lowercase();
        }
        let h = shared_label(r, "s", i, 0.85, 0.3, 9);
        seqs.push((h, s));
    }
    if r.chance(0.03) {
        let (a, b) = big(r);
        seqs.push(("big0".into(), a));
        seqs.push(("big1".into(), b));
    }
    w(t, "in.fa", &fasta_all(&seqs))?;
    let mut o: Vec<String> = Vec::new();
    if r.chance(0.3) {
        o.push("acceptall".into());
    } else {
        o.extend([
            "id".into(),
            r.choice(&["0", "0.5", "0.7", "0.8", "0.9", "0.97", "1"])
                .into(),
        ]);
    }
    let mut table: Vec<(&str, &[&str])> = vec![
        ("iddef", &["0", "1", "3", "4"]),
        ("weak_id", &["0.3", "0.6"]),
        ("maxhits", &["1", "2", "4"]),
    ];
    table.extend_from_slice(CRITERIA);
    options(r, &mut o, 0.07, &table);
    flags(
        r,
        &mut o,
        0.07,
        &[
            "leftjust",
            "rightjust",
            "self",
            "selfid",
            "top_hits_only",
            "n_mismatch",
        ],
    );
    scores(r, &mut o, 0.2);
    gaps(
        r,
        &mut o,
        0.25,
        &[
            "10I/1E",
            "20I/2E",
            "5Q/8T",
            "*E/20I",
            "3I/1LQ/2RT",
            "0E/20I",
            "1I/0E",
            "*I/2E",
        ],
        GAPEXT,
    );
    params(t, &o)
}

pub fn dbsrch(t: &Path, r: &mut Rng) -> io::Result<()> {
    // queries and a database of sequences related by substitutions, indels,
    // overhangs and repeats, in mixed case with IUPAC codes, empty, short or
    // random ones, now and then long ones for the fallback aligner;
    // abundances and labels shared between queries and targets; and
    // options: identity, limits, words, strands, criteria, scores, gap
    // penalties; a UDB database now and then (made by mkudb in the case)
    let len: usize = r.choice(&[20, 60, 150, 400]);
    let n = r.randint(1, 3) as usize;
    let tem = templates(r, len, n);
    let one = |r: &mut Rng| -> String {
        let k = r.random();
        let mut s = if k < 0.04 {
            String::new()
        } else if k < 0.1 {
            let n = r.randint(1, len as i64) as usize;
            rnd(r, n)
        } else {
            piece(r, &tem, 4, &[0.0, 0.01, 0.05, 0.15, 0.3], 0.2, 20)
        };
        if r.chance(0.15) {
            s = rc(&s);
        }
        if r.chance(0.1) {
            s = sprinkle(r, &s, 0.05);
        }
        if r.chance(0.1) {
            s = s.to_lowercase();
        } else if r.chance(0.1) {
            s = lower_span(r, &s);
        }
        s
    };
    let mut db: Vec<(String, String)> = (0..r.choice(&[0, 1, 3, 10, 30, 60, 120]))
        .map(|i| {
            let h = shared_label(r, "t", i, 0.9, 0.3, 9);
            (h, one(r))
        })
        .collect();
    let mut qs: Vec<(String, String)> = (0..r.randint(1, 15) as usize)
        .map(|i| {
            let h = shared_label(r, "q", i, 0.9, 0.3, 9);
            (h, one(r))
        })
        .collect();
    if r.chance(0.03) {
        let (a, b) = big(r);
        db.push(("big0".into(), a));
        qs.push(("bigq".into(), b));
    }
    w(t, "db.fa", &fasta_all(&db))?;
    w(t, "q.fa", &fasta_all(&qs))?;
    let mut o: Vec<String> = vec![
        "id".into(),
        r.choice(&["0", "0.5", "0.7", "0.8", "0.9", "0.97", "1"])
            .into(),
    ];
    options(
        r,
        &mut o,
        0.3,
        &[
            ("maxaccepts", &["0", "1", "2", "5", "20"]),
            ("maxrejects", &["0", "1", "3", "8", "32", "100"]),
            ("wordlength", &["3", "4", "6", "8", "10", "12", "15"]),
            ("minwordmatches", &["0", "1", "3", "10", "30"]),
            ("strand", &["both"]),
            ("iddef", &["0", "1", "3", "4"]),
            ("weak_id", &["0.3", "0.6", "0.85"]),
            ("maxhits", &["1", "2", "4"]),
        ],
    );
    options(r, &mut o, 0.05, CRITERIA);
    flags(
        r,
        &mut o,
        0.07,
        &[
            "leftjust",
            "rightjust",
            "self",
            "selfid",
            "top_hits_only",
            "n_mismatch",
        ],
    );
    scores(r, &mut o, 0.2);
    gaps(r, &mut o, 0.2, GAPOPEN, GAPEXT);
    if !db.is_empty() && r.chance(0.15) {
        let k = r.choice(&["4", "6", "8", "11"]);
        o.splice(0..0, ["udb".to_string(), k.to_string()]);
    }
    params(t, &o)
}

pub fn clust(t: &Path, r: &mut Rng) -> io::Result<()> {
    // sequences related by substitutions, indels, overhangs and repeats to
    // a few templates, reverse complements, mixed case, IUPAC codes, empty
    // and random ones, now and then long ones for the fallback aligner;
    // abundances and shared labels; options: identity, limits, words,
    // strands, size order, abundances in and out, criteria, scores, gap
    // penalties; and the mode: input order (smallmem), or sorted by length
    // (fast) or size
    let len: usize = r.choice(&[20, 60, 150, 400]);
    let n = r.randint(1, 5) as usize;
    let tem = templates(r, len, n);
    let mut seqs = Vec::new();
    for i in 0..r.choice(&[0, 1, 5, 20, 60, 150]) {
        let k = r.random();
        let mut s = if k < 0.03 {
            String::new()
        } else if k < 0.08 {
            let n = r.randint(1, len as i64) as usize;
            rnd(r, n)
        } else {
            piece(r, &tem, 5, &[0.0, 0.01, 0.03, 0.06, 0.15], 0.1, 15)
        };
        if r.chance(0.1) {
            s = rc(&s);
        }
        if r.chance(0.05) {
            s = sprinkle(r, &s, 0.05);
        }
        if r.chance(0.05) {
            s = s.to_lowercase();
        } else if r.chance(0.08) {
            s = lower_span(r, &s);
        }
        let h = shared_label(r, "s", i, 0.92, 0.5, 30);
        seqs.push((h, s));
    }
    if r.chance(0.03) {
        let (a, b) = big(r);
        seqs.push(("big0".into(), a));
        seqs.push(("big1".into(), b));
    }
    w(t, "in.fa", &fasta_all(&seqs))?;
    let mode = r.choice(&["smallmem", "smallmem", "fast", "size"]);
    let id = r.choice(&["0.5", "0.7", "0.8", "0.9", "0.95", "0.97", "0.99", "1"]);
    let mut o: Vec<String> = vec![mode.into(), "id".into(), id.into()];
    options(
        r,
        &mut o,
        0.25,
        &[
            ("maxaccepts", &["0", "1", "2", "5", "20"]),
            ("maxrejects", &["0", "1", "3", "8", "32", "100"]),
            ("wordlength", &["3", "4", "6", "8", "10", "12", "15"]),
            ("minwordmatches", &["0", "1", "3", "10", "30"]),
            ("strand", &["both"]),
            ("iddef", &["0", "1", "3", "4"]),
            ("weak_id", &["0.3", "0.6"]),
        ],
    );
    flags(r, &mut o, 0.15, &["sizeorder", "nosizein", "nosizeout"]);
    options(r, &mut o, 0.04, CRITERIA);
    flags(
        r,
        &mut o,
        0.05,
        &["leftjust", "rightjust", "self", "selfid", "n_mismatch"],
    );
    scores(r, &mut o, 0.15);
    gaps(r, &mut o, 0.15, GAPOPEN, GAPEXT);
    params(t, &o)
}

pub fn unoise(t: &Path, r: &mut Rng) -> io::Result<()> {
    // abundant templates and their variants by substitutions and indels, of
    // every abundance, reverse complements, mixed case, IUPAC codes, empty
    // and random ones, shared abundances (ties); options: --minsize (sqfilt
    // before), alpha, identity, limits, words, strands, size order,
    // abundances in and out, criteria, scores
    let len: usize = r.choice(&[30, 80, 150, 300]);
    let tem: Vec<String> = (0..r.randint(1, 6)).map(|_| rnd(r, len)).collect();
    let mut seqs: Vec<(String, String)> = tem
        .iter()
        .enumerate()
        .map(|(i, s)| {
            (
                format!("t{i};size={}", r.choice(&[50, 100, 300, 1000, 5000])),
                s.clone(),
            )
        })
        .collect();
    for i in 0..r.choice(&[0, 5, 20, 60, 150]) {
        let k = r.random();
        let mut s = if k < 0.03 {
            String::new()
        } else if k < 0.07 {
            let n = r.randint(1, len as i64) as usize;
            rnd(r, n)
        } else {
            let tm = r.choice(&tem);
            let n = r.choice(&[0, 1, 1, 1, 2, 2, 3, 4, 6, 10]);
            mutate(r, &tm, n)
        };
        if r.chance(0.08) {
            s = rc(&s);
        }
        if r.chance(0.04) {
            s = sprinkle(r, &s, 0.05);
        }
        if r.chance(0.04) {
            s = s.to_lowercase();
        } else if r.chance(0.06) {
            s = lower_span(r, &s);
        }
        let h = if r.chance(0.95) {
            format!("s{i}")
        } else {
            format!("dup{}", r.randrange(2))
        };
        let size = r.choice(&[1, 1, 2, 3, 4, 6, 8, 10, 15, 25, 40, 80, 200]);
        seqs.push((format!("{h};size={size}"), s));
    }
    r.shuffle(&mut seqs);
    w(t, "in.fa", &fasta_all(&seqs))?;
    let mut o: Vec<String> = vec![
        "minsize".into(),
        r.choice(&["1", "1", "2", "4", "8"]).into(),
    ];
    options(
        r,
        &mut o,
        0.25,
        &[
            ("unoise_alpha", &["0", "0.5", "1", "2", "3", "5"]),
            ("weak_id", &["0.5", "0.8", "0.9", "0.95", "0.99"]),
            ("maxaccepts", &["0", "1", "2", "5"]),
            ("maxrejects", &["0", "1", "3", "8", "32"]),
            ("wordlength", &["3", "5", "8", "12"]),
            ("minwordmatches", &["0", "1", "5", "20"]),
            ("strand", &["both"]),
            ("iddef", &["0", "1", "3", "4"]),
        ],
    );
    flags(r, &mut o, 0.15, &["sizeorder", "nosizein", "nosizeout"]);
    options(
        r,
        &mut o,
        0.05,
        &[
            ("maxdiffs", &["0", "2", "10"]),
            ("maxgaps", &["0", "1"]),
            ("query_cov", &["0.5", "0.9"]),
            ("target_cov", &["0.9"]),
            ("maxid", &["0.99"]),
            ("mincols", &["20"]),
        ],
    );
    flags(
        r,
        &mut o,
        0.04,
        &["leftjust", "rightjust", "self", "selfid", "n_mismatch"],
    );
    if r.chance(0.12) {
        let (m, x) = (r.choice(&["1", "3"]), r.choice(&["-1", "-2", "-6"]));
        o.extend(["match".into(), m.into(), "mismatch".into(), x.into()]);
    }
    gaps(
        r,
        &mut o,
        0.1,
        &["10I/1E", "5Q/8T", "*I/2E"],
        &["1I/1E", "2I/3E"],
    );
    params(t, &o)
}

/// `s` with `k` substitutions at distinct places, each to another base.
fn subst_distinct(r: &mut Rng, s: &str, k: usize) -> String {
    let mut b = s.as_bytes().to_vec();
    for p in r.sample(b.len(), k) {
        let others: Vec<u8> = b"ACGT".iter().copied().filter(|&c| c != b[p]).collect();
        b[p] = r.choice(&others);
    }
    String::from_utf8(b).unwrap()
}

/// The finishing touches of a chimera case's sequence: more edits, IUPAC
/// codes, lower case.
fn touch(r: &mut Rng, mut s: String, edits: f64, edit_counts: &[usize], iupac: f64) -> String {
    if r.chance(edits) {
        let n = r.choice(edit_counts);
        s = mutate(r, &s, n);
    }
    if r.chance(iupac) {
        s = sprinkle(r, &s, 0.03);
    }
    if r.chance(0.05) {
        s = s.to_lowercase();
    } else if r.chance(0.08) {
        s = lower_span(r, &s);
    }
    s
}

pub fn uchime(t: &Path, r: &mut Rng) -> io::Result<()> {
    // parents, some a few substitutions from another, and chimeras of two or
    // three of them, with substitutions and indels, in mixed case, with
    // IUPAC codes, empty, short and random sequences, abundances (parents
    // the more abundant, now and then not); and options: the command
    // (uchime, uchime2, uchime3 de novo, or against a database), scores and
    // thresholds, abundance skew, alignment width
    let len: usize = r.choice(&[40, 100, 200, 400]);
    let mut par: Vec<String> = (0..r.randint(2, 6)).map(|_| rnd(r, len)).collect();
    if r.chance(0.5) {
        let core = rnd(r, len / 3);
        par = par
            .into_iter()
            .map(|p| {
                if r.chance(0.5) {
                    format!("{}{core}{}", &p[..len / 3], &p[2 * len / 3..])
                } else {
                    p
                }
            })
            .collect();
    }
    for _ in 0..r.randint(0, 3) {
        let p = r.choice(&par);
        let k = r.choice(&[2, 5, 10, 20]);
        par.push(subst_distinct(r, &p, k));
    }
    let mut seqs: Vec<(String, String)> = par
        .iter()
        .enumerate()
        .map(|(i, p)| {
            (
                format!("p{i};size={}", r.choice(&[100, 200, 500, 1000])),
                p.clone(),
            )
        })
        .collect();
    for i in 0..r.choice(&[0, 3, 10, 30, 60]) {
        let k = r.random();
        let s = if k < 0.6 {
            let ab = r.sample(par.len(), 2);
            let c = r.randint(1, len as i64 - 1) as usize;
            format!("{}{}", &par[ab[0]][..c], &par[ab[1]][c..])
        } else if k < 0.75 && par.len() > 2 {
            let abd = r.sample(par.len(), 3);
            let mut cs: Vec<usize> = r.sample(len - 1, 2).into_iter().map(|c| c + 1).collect();
            cs.sort();
            format!(
                "{}{}{}",
                &par[abd[0]][..cs[0]],
                &par[abd[1]][cs[0]..cs[1]],
                &par[abd[2]][cs[1]..]
            )
        } else if k < 0.9 {
            let p = r.choice(&par);
            let n = r.choice(&[0, 1, 3, 8]);
            mutate(r, &p, n)
        } else if k < 0.95 {
            let n = r.randint(0, len as i64) as usize;
            rnd(r, n)
        } else {
            let n = r.randint(0, 5) as usize;
            rnd(r, n)
        };
        let s = touch(r, s, 0.25, &[1, 2, 4], 0.05);
        let h = if r.chance(0.95) {
            format!("c{i}")
        } else {
            format!("dup{}", r.randrange(2))
        };
        let size = r.choice(&[1, 2, 3, 5, 10, 50, 300, 2000]);
        seqs.push((format!("{h};size={size}"), s));
    }
    r.shuffle(&mut seqs);
    let mode = r.choice(&[
        "uchime_denovo",
        "uchime_denovo",
        "uchime2_denovo",
        "uchime3_denovo",
        "uchime_ref",
    ]);
    if mode == "uchime_ref" {
        let mut refs: Vec<(String, String)> = seqs
            .iter()
            .filter(|x| x.0.starts_with('p'))
            .cloned()
            .collect();
        if r.chance(0.3) {
            refs.push(("rX".into(), rnd(r, len)));
        }
        let qs: Vec<(String, String)> = seqs
            .iter()
            .filter(|x| !x.0.starts_with('p') || r.chance(0.2))
            .cloned()
            .collect();
        w(t, "db.fa", &fasta_all(&refs))?;
        w(t, "in.fa", &fasta_all(&qs))?;
    } else {
        w(t, "in.fa", &fasta_all(&seqs))?;
    }
    let mut o: Vec<String> = vec![mode.into()];
    let mut table: Vec<(&str, &[&str])> = vec![
        ("xn", &["1.5", "4", "16"]),
        ("dn", &["0.5", "3"]),
        ("alignwidth", &["0", "30", "60"]),
    ];
    if mode == "uchime_denovo" || mode == "uchime_ref" {
        table.extend([
            ("minh", &["0.05", "0.1", "0.5", "1"][..]),
            ("mindiv", &["0.2", "2", "10"]),
            ("mindiffs", &["1", "5", "10"]),
        ]);
    }
    if mode != "uchime_ref" {
        table.push(("abskew", &["1", "1.5", "4", "16"]));
    }
    options(r, &mut o, 0.25, &table);
    params(t, &o)
}

pub fn mchime(t: &Path, r: &mut Rng) -> io::Result<()> {
    // long parents, some a few substitutions from another, and chimeras of
    // two to five of them, exact or with substitutions and indels, in mixed
    // case, with IUPAC codes, empty, short and random sequences,
    // abundances; and options: parts, parents, segment length, differences,
    // abundance skew, alignment width
    let len: usize = r.choice(&[60, 200, 500, 1000]);
    let mut par: Vec<String> = (0..r.randint(2, 7)).map(|_| rnd(r, len)).collect();
    for _ in 0..r.randint(0, 3) {
        let p = r.choice(&par);
        let k = r.choice(&[3, 10, 30]);
        par.push(subst_distinct(r, &p, k));
    }
    let mut seqs: Vec<(String, String)> = par
        .iter()
        .enumerate()
        .map(|(i, p)| {
            (
                format!("p{i};size={}", r.choice(&[100, 200, 500, 1000])),
                p.clone(),
            )
        })
        .collect();
    for i in 0..r.choice(&[0, 3, 10, 30]) {
        let k = r.random();
        let s = if k < 0.7 {
            let n = par.len().min(r.choice(&[2, 2, 3, 3, 4, 5]));
            let ps = r.sample(par.len(), n);
            let mut b: Vec<usize> = r
                .sample(len - 1, n - 1)
                .into_iter()
                .map(|c| c + 1)
                .collect();
            b.sort();
            b.insert(0, 0);
            b.push(len);
            (0..n).map(|j| &par[ps[j]][b[j]..b[j + 1]]).collect()
        } else if k < 0.85 {
            let p = r.choice(&par);
            let n = r.choice(&[0, 1, 3]);
            mutate(r, &p, n)
        } else if k < 0.95 {
            let n = r.randint(0, len as i64) as usize;
            rnd(r, n)
        } else {
            let n = r.randint(0, 5) as usize;
            rnd(r, n)
        };
        let s = touch(r, s, 0.2, &[1, 2], 0.04);
        let h = if r.chance(0.95) {
            format!("c{i}")
        } else {
            format!("dup{}", r.randrange(2))
        };
        let size = r.choice(&[1, 2, 3, 5, 10, 50, 300, 2000]);
        seqs.push((format!("{h};size={size}"), s));
    }
    r.shuffle(&mut seqs);
    w(t, "in.fa", &fasta_all(&seqs))?;
    let mut o = Vec::new();
    options(
        r,
        &mut o,
        0.25,
        &[
            ("chimeras_parts", &["2", "3", "5", "10", "100"]),
            ("chimeras_parents_max", &["2", "3", "4", "6", "20"]),
            ("chimeras_length_min", &["1", "5", "20", "100"]),
            ("chimeras_diff_pct", &["0.5", "2", "10", "50"]),
            ("abskew", &["1", "2", "8"]),
            ("alignwidth", &["0", "40", "80"]),
        ],
    );
    params(t, &o)
}
