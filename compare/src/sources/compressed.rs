//! Compressed input. The generators' files are gzipped or bzip2'd. vsearch
//! reads them itself (it knows them by their magic numbers), and the
//! translated pipeline decompresses them with gzip -dc or bzip2 -dc.
//! Malformed records inside the compressed files, and files cut short, are
//! compared too.

use std::io;
use std::path::Path;
use std::process::Command;

use super::{Env, Source, argv, generate, p, seeds};
use crate::cases::Case;
use crate::rng::Rng;

pub fn sources() -> Vec<Source> {
    vec![Source {
        name: "compressed",
        seeds: || seeds(1, 100),
        about: "gzip and bzip2 input: streaming and whole-file tools, pairs, databases, malformed records, \
                files cut short",
        make: compressed,
    }]
}

/// The file compressed by `tool` next to it, `.gz` or `.bz2`.
fn pack(path: &Path, tool: &str) -> io::Result<String> {
    let ext = if tool == "gzip" { "gz" } else { "bz2" };
    let out = format!("{}.{ext}", p(path));
    let data = Command::new(tool).arg("-c").arg(path).output()?;
    if !data.status.success() {
        return Err(io::Error::other(format!("{tool} -c {} failed", p(path))));
    }
    std::fs::write(&out, data.stdout)?;
    Ok(out)
}

fn compressed(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    let mut r = Rng::new(seed);
    let mut kind = || if r.chance(0.5) { "gzip" } else { "bzip2" };
    for g in ["derep", "mergpa", "dbsrch", "fuzz"] {
        let d = dir.join(g);
        std::fs::create_dir_all(&d)?;
        generate(g, &d, seed, &[])?;
    }
    let (d, m, s, f) = (
        dir.join("derep"),
        dir.join("mergpa"),
        dir.join("dbsrch"),
        dir.join("fuzz"),
    );
    let q = pack(&d.join("q.fq"), kind())?;
    let a = pack(&d.join("a.fa"), kind())?;
    let mut out = vec![
        Case::new(
            "uniques fq",
            argv(&format!("--fastx_uniques {q} --quiet --sizein --sizeout --minseqlength 0 --notrunclabels --fastqout {{O}}/o")),
        ),
        Case::new(
            "filter fq",
            argv(&format!(
                "--fastx_filter {q} --quiet --notrunclabels --fastq_stripleft 2 --fastq_maxee 1 --fastqout {{O}}/k \
                 --fastqout_discarded {{O}}/d"
            )),
        )
        .success_outputs_only(),
        Case::new(
            "revcomp fa",
            argv(&format!("--fastx_revcomp {a} --quiet --notrunclabels --fasta_width 0 --fastaout {{O}}/o")),
        ),
        Case::new(
            "sortbylength fa",
            argv(&format!(
                "--sortbylength {a} --quiet --notrunclabels --minseqlength 0 --fasta_width 0 --output {{O}}/o"
            )),
        ),
        Case::new("fasta at 80", argv(&format!("--fastx_revcomp {q} --quiet --notrunclabels --fastaout {{O}}/o"))),
        Case::new("chars", argv(&format!("--fastq_chars {q} --quiet --log {{O}}/log"))),
    ];
    let (fw, rv) = (
        pack(&m.join("f.fq"), kind())?,
        pack(&m.join("r.fq"), kind())?,
    );
    out.push(
        Case::new(
            "mergepairs",
            argv(&format!(
                "--fastq_mergepairs {fw} --reverse {rv} --quiet --fastqout {{O}}/m --fastqout_notmerged_fwd {{O}}/n \
                 --eetabbedout {{O}}/e --log {{O}}/log"
            )),
        )
        .errors("mergpa")
        .success_outputs_only(),
    );
    let (sq, sdb) = (
        pack(&s.join("q.fa"), kind())?,
        pack(&s.join("db.fa"), kind())?,
    );
    out.push(
        Case::new(
            "search",
            argv(&format!(
                "--usearch_global {sq} --db {sdb} --id 0.7 --uc {{O}}/uc --uc_allhits --quiet --notrunclabels \
                 --qmask soft --dbmask soft --minseqlength 0 --maxseqlength 100000 --threads 1"
            )),
        )
        .success_only(&["uc"]),
    );
    // malformed records inside compressed files
    let (ff, fa) = (
        pack(&f.join("f.fq"), kind())?,
        pack(&f.join("f.fa"), kind())?,
    );
    for (name, v) in [
        (
            "fuzz filter fq",
            format!("--fastx_filter {ff} --quiet --notrunclabels --fastqout {{O}}/o"),
        ),
        (
            "fuzz revcomp fa",
            format!(
                "--fastx_revcomp {fa} --quiet --notrunclabels --fasta_width 0 --fastaout {{O}}/o"
            ),
        ),
        (
            "fuzz uniques fa",
            format!(
                "--fastx_uniques {fa} --quiet --notrunclabels --minseqlength 0 --sizein --sizeout --fasta_width 0 \
                 --fastaout {{O}}/o"
            ),
        ),
    ] {
        out.push(
            Case::new(name, argv(&v))
                .errors("exact")
                .warn(&["stripped"]),
        );
    }
    // a compressed file cut short: the exit code only, as what comes out
    // before the error is the decompressor's (macOS gzip -dc writes
    // nothing, bzip2 -dc what it decoded)
    let data = std::fs::read(&q)?;
    let ext = if q.ends_with(".gz") { "gz" } else { "bz2" };
    let cut = format!("{q}.cut.{ext}");
    let n = r.randint(1, (data.len() as i64 - 1).max(1)) as usize;
    std::fs::write(&cut, &data[..n])?;
    out.push(
        Case::new(
            "cut short",
            argv(&format!(
                "--fastx_revcomp {cut} --quiet --notrunclabels --fastqout {{O}}/o"
            )),
        )
        .success_outputs_only(),
    );
    Ok(out)
}
