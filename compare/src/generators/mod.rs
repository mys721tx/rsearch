//! The input generators. Each writes one case's input files (and a `params`
//! file of drawn options where the case needs one) into a directory, from a
//! seed. A seed always writes the same files.

use std::fs;
use std::io;
use std::path::Path;

use crate::rng::Rng;

mod ana;
mod fuzz;
mod pre;

/// The generators, by name.
pub const NAMES: &[&str] = &[
    "allpr", "clust", "conv", "dbornt", "dbsrch", "derep", "exsrch", "fqchar", "fqstat", "fuzz",
    "getseq", "hitfmt", "joinpa", "lcatax", "lcmask", "mchime", "mergpa", "mkudb", "otutab",
    "redig", "relab", "rerep", "revcomp", "scrmbl", "sff", "shufl", "sintax", "sortby", "sqfilt",
    "sqtrim", "subseq", "subsmp", "syncpa", "uchime", "udbinfo", "udbstat", "unoise",
];

/// Writes generator `name`'s inputs for `seed` into `dir`. `mkudb` takes
/// an optional word length, `udbstat` the path of the mkudb tool.
pub fn write(name: &str, dir: &Path, seed: u64, extra: &[&str]) -> io::Result<()> {
    let r = &mut Rng::new(seed);
    match name {
        "allpr" => ana::allpr(dir, r),
        "clust" => ana::clust(dir, r),
        "conv" => pre::conv(dir, r),
        "dbornt" => ana::dbornt(dir, r),
        "dbsrch" => ana::dbsrch(dir, r),
        "derep" => pre::derep(dir, r),
        "exsrch" => ana::exsrch(dir, r),
        "fqchar" => pre::fqchar(dir, r),
        "fqstat" => pre::fqstat(dir, r),
        "fuzz" => fuzz::fuzz(dir, r),
        "getseq" => pre::getseq(dir, r),
        "hitfmt" => ana::hitfmt(dir, r),
        "joinpa" => pre::joinpa(dir, r),
        "lcatax" => ana::lcatax(dir, r),
        "lcmask" => pre::lcmask(dir, r),
        "mchime" => ana::mchime(dir, r),
        "mergpa" => pre::mergpa(dir, r),
        "mkudb" => {
            let k = extra
                .first()
                .map(|k| k.parse::<u32>().map_err(io::Error::other))
                .transpose()?;
            ana::mkudb(dir, r, k)
        }
        "otutab" => ana::otutab(dir, r),
        "redig" => pre::redig(dir, r),
        "relab" => pre::relab(dir, r),
        "rerep" => pre::rerep(dir, r),
        "revcomp" => pre::revcomp(dir, r),
        "scrmbl" => pre::scrmbl(dir, r),
        "sff" => pre::sff(dir, r),
        "shufl" => pre::shufl(dir, r),
        "sintax" => ana::sintax(dir, r),
        "sortby" => pre::sortby(dir, r),
        "sqfilt" => pre::sqfilt(dir, r),
        "sqtrim" => pre::sqtrim(dir, r),
        "subseq" => pre::subseq(dir, r),
        "subsmp" => pre::subsmp(dir, r),
        "syncpa" => pre::syncpa(dir, r),
        "uchime" => ana::uchime(dir, r),
        "udbinfo" => ana::udbinfo(dir, r),
        "udbstat" => {
            let mkudb = extra
                .first()
                .ok_or_else(|| io::Error::other("udbstat needs the mkudb tool"))?;
            ana::udbstat(dir, r, Path::new(mkudb))
        }
        "unoise" => ana::unoise(dir, r),
        _ => Err(io::Error::other(format!("no generator {name}"))),
    }
}

/// Text built record by record.
#[derive(Default)]
pub struct Text(pub String);

impl Text {
    pub fn new() -> Text {
        Text(String::new())
    }

    pub fn push(&mut self, s: &str) {
        self.0.push_str(s);
    }
}

pub fn fasta(h: &str, s: &str) -> String {
    format!(">{h}\n{s}\n")
}

pub fn fastq(h: &str, s: &str, q: &str) -> String {
    format!("@{h}\n{s}\n+\n{q}\n")
}

/// Records as FASTA.
pub fn fasta_all(recs: &[(String, String)]) -> String {
    recs.iter().map(|(h, s)| fasta(h, s)).collect()
}

/// Records as FASTQ, every score 40.
pub fn fastq_i(recs: &[(String, String)]) -> String {
    recs.iter()
        .map(|(h, s)| fastq(h, s, &"I".repeat(s.len())))
        .collect()
}

/// The reverse complement, IUPAC codes and case kept, U read as T.
pub fn rc(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|c| {
            let x = match c.to_ascii_uppercase() {
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
                other => other,
            };
            (if c.is_ascii_lowercase() {
                x.to_ascii_lowercase()
            } else {
                x
            }) as char
        })
        .collect()
}

pub fn w(dir: &Path, name: &str, text: &str) -> io::Result<()> {
    fs::write(dir.join(name), text)
}

pub fn wb(dir: &Path, name: &str, data: &[u8]) -> io::Result<()> {
    fs::write(dir.join(name), data)
}

/// `n` random bases.
pub fn rnd(r: &mut Rng, n: usize) -> String {
    r.chars("ACGT", n)
}

/// `s[a..a + n]`, cut at its end.
pub fn slice(s: &str, a: usize, n: usize) -> &str {
    let a = a.min(s.len());
    &s[a..(a + n).min(s.len())]
}

/// `s` with a random stretch in lower case.
pub fn lower_span(r: &mut Rng, s: &str) -> String {
    let a = r.randrange(s.len() + 1);
    let b = a + r.randrange(s.len() - a + 1);
    format!("{}{}{}", &s[..a], s[a..b].to_lowercase(), &s[b..])
}

/// `s` with each residue an ambiguous code with probability `p`.
pub fn sprinkle(r: &mut Rng, s: &str, p: f64) -> String {
    s.chars()
        .map(|c| {
            if r.chance(p) {
                r.ch("NRYKM") as char
            } else {
                c
            }
        })
        .collect()
}

/// `s` with `k` substitutions at random places (some may repeat).
pub fn substitute(r: &mut Rng, s: &str, k: usize) -> String {
    let mut s = s.as_bytes().to_vec();
    for _ in 0..k {
        if !s.is_empty() {
            let p = r.randrange(s.len());
            s[p] = r.ch("ACGT");
        }
    }
    String::from_utf8(s).unwrap()
}

/// `s` with `k` edits: substitutions mostly, deletions and insertions.
pub fn mutate(r: &mut Rng, s: &str, k: usize) -> String {
    let mut s = s.as_bytes().to_vec();
    for _ in 0..k {
        if s.is_empty() {
            break;
        }
        let p = r.randrange(s.len());
        let x = r.random();
        if x < 0.7 {
            s[p] = r.ch("ACGT");
        } else if x < 0.85 {
            s.remove(p);
        } else {
            let c = r.ch("ACGT");
            s.insert(p, c);
        }
    }
    String::from_utf8(s).unwrap()
}

/// `s` with each residue deleted, substituted or followed by an insertion
/// at a third of `rate` each (and, with `iupac`, made ambiguous at a tenth).
pub fn edit(r: &mut Rng, s: &str, rate: f64, iupac: bool) -> String {
    let mut out = String::new();
    for c in s.chars() {
        let k = r.random();
        if k < rate / 3.0 {
            continue;
        }
        if k < 2.0 * rate / 3.0 {
            out.push(r.ch("ACGT") as char);
        } else if k < rate {
            out.push(c);
            let n = r.randint(1, if iupac { 4 } else { 6 }) as usize;
            out += &rnd(r, n);
        } else if iupac && k < rate * 1.1 {
            out.push(r.ch("NRYKM") as char);
        } else {
            out.push(c);
        }
    }
    out
}

/// Adds each named option with probability `p`, with a random value.
pub fn options(r: &mut Rng, o: &mut Vec<String>, p: f64, table: &[(&str, &[&str])]) {
    for (name, vals) in table {
        if r.chance(p) {
            o.push(name.to_string());
            o.push(r.choice(vals).to_string());
        }
    }
}

/// Adds each named flag with probability `p`.
pub fn flags(r: &mut Rng, o: &mut Vec<String>, p: f64, names: &[&str]) {
    for name in names {
        if r.chance(p) {
            o.push(name.to_string());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn helpers() {
        assert_eq!(rc("ACGTUacgun"), "nacgtAACGT");
        assert_eq!(rc("RYKMBVDHSW"), "WSDHBVKMRY");
        assert_eq!(slice("ACGT", 2, 10), "GT");
        assert_eq!(slice("ACGT", 9, 1), "");
        let r = &mut Rng::new(1);
        for _ in 0..50 {
            let s = lower_span(r, "ACGTACGT");
            assert_eq!(s.to_uppercase(), "ACGTACGT");
            assert_eq!(mutate(r, "", 3), "");
            assert_eq!(edit(r, "ACGT", 0.0, true), "ACGT");
        }
    }

    /// Every generator writes something for a few seeds, and the same
    /// files for the same seed.
    #[test]
    fn every_generator() {
        let root = std::env::temp_dir().join(format!("rscompare-gen-{}", std::process::id()));
        for name in NAMES {
            if *name == "udbstat" {
                continue;
            }
            for seed in 0..8 {
                let a = root.join(format!("{name}-{seed}-a"));
                let b = root.join(format!("{name}-{seed}-b"));
                for d in [&a, &b] {
                    fs::create_dir_all(d).unwrap();
                    write(name, d, seed, &[]).unwrap();
                }
                let mut files: Vec<_> = fs::read_dir(&a)
                    .unwrap()
                    .map(|e| e.unwrap().file_name())
                    .collect();
                files.sort();
                assert!(!files.is_empty(), "{name} {seed}");
                for f in files {
                    assert_eq!(
                        fs::read(a.join(&f)).unwrap(),
                        fs::read(b.join(&f)).unwrap(),
                        "{name} {seed} {f:?}"
                    );
                }
            }
        }
        fs::remove_dir_all(&root).unwrap();
    }
}
