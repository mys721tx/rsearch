//! Malformed FASTA and FASTQ: odd bytes in sequences, headers and scores,
//! CR LF line ends, wrapped lines, missing final newlines, cut-short files.

use std::io;
use std::path::Path;

use super::wb;
use crate::rng::Rng;

const ODD: &[&[u8]] = &[
    b"1", b"*", b"-", b".", b" ", b"\t", b"\x0b", b"\x0c", b"\r", b"\x01", b"\x7f", b"\xe9", b"x",
    b"U", b"u",
];

fn seq(r: &mut Rng, n: usize, rare: f64) -> Vec<u8> {
    let mut out = Vec::new();
    for _ in 0..n {
        if r.chance(rare) {
            out.extend(r.choice(ODD));
        } else {
            out.push(r.ch("ACGTNacgtRYn"));
        }
    }
    out
}

fn header(r: &mut Rng, i: usize) -> Vec<u8> {
    let mut h = format!("r{i}").into_bytes();
    let k = r.random();
    if k < 0.1 {
        h.extend(b" desc x");
    } else if k < 0.15 {
        h.extend(b"\tt");
    } else if k < 0.17 {
        h.extend(b"\rcut");
    } else if k < 0.18 {
        h.extend(r.choice(&[b"\x0b", b"\x7f", b"\x01"]));
        h.push(b'z');
    } else if k < 0.2 {
        h.clear();
    }
    h
}

/// `s` with a line break somewhere inside, a fifth of the time.
fn wrap(r: &mut Rng, s: &[u8]) -> Vec<u8> {
    if s.len() > 3 && r.chance(0.2) {
        let k = r.randint(1, s.len() as i64 - 1) as usize;
        return [&s[..k], b"\n", &s[k..]].concat();
    }
    s.to_vec()
}

fn strip_newlines(s: &mut Vec<u8>) {
    while s.last() == Some(&b'\n') {
        s.pop();
    }
}

pub fn fuzz(t: &Path, r: &mut Rng) -> io::Result<()> {
    let rare = r.choice(&[0.0, 0.002, 0.02]);
    let (mut fa, mut fq) = (Vec::new(), Vec::new());
    for i in 0..r.randint(1, 15) as usize {
        let n = r.randint(0, 20) as usize;
        let s = seq(r, n, rare);
        let h = header(r, i);
        let eol: &[u8] = if r.chance(0.1) { b"\r\n" } else { b"\n" };
        fa.push(b'>');
        fa.extend(&h);
        fa.extend(eol);
        fa.extend(wrap(r, &s));
        fa.extend(eol);
        // quality: mostly printable, rarely bad
        let q: Vec<u8> = (0..s.len())
            .map(|_| {
                if r.random() >= rare {
                    r.randint(33, 73) as u8
                } else {
                    r.choice(&[32, 9, 127, 13])
                }
            })
            .collect();
        let plus = r.weighted(
            &[b"".to_vec(), h.clone(), b"other".to_vec()],
            &[0.85, 0.1, 0.05],
        );
        fq.push(b'@');
        fq.extend(&h);
        fq.extend(eol);
        fq.extend(wrap(r, &s));
        fq.extend(eol);
        fq.push(b'+');
        fq.extend(plus);
        fq.extend(eol);
        fq.extend(wrap(r, &q));
        fq.extend(eol);
    }
    if r.chance(0.15) {
        strip_newlines(&mut fa);
    }
    if r.chance(0.15) {
        strip_newlines(&mut fq);
    }
    if r.chance(0.05) {
        fa.extend(b">end");
    }
    if r.chance(0.05) {
        fq.truncate(r.randint(1, fq.len() as i64) as usize);
    }
    wb(t, "f.fa", &fa)?;
    wb(t, "f.fq", &fq)
}
