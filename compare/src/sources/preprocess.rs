//! Cases for the preprocessing tools, from the former zsh runners with
//! the same generators, option sets and seed ranges.

use std::io;
use std::path::Path;

use super::{Env, Source, argv, generate, p, read, seeds};
use crate::cases::Case;
use crate::rng::Rng;
use crate::translate::Step;

pub fn sources() -> Vec<Source> {
    vec![
        Source {
            name: "derep",
            seeds: || seeds(1, 5),
            about: "qualities, encodings, uc, tabbedout, --same-label, --prefix",
            make: derep,
        },
        Source {
            name: "conv",
            seeds: || seeds(1, 3),
            about: "FASTQ and FASTA: encodings, Solexa, ranges, labels, errors",
            make: conv,
        },
        Source {
            name: "sff",
            seeds: || seeds(1, 60),
            about: "SFF: clipping, encodings, labels, truncated files",
            make: sff,
        },
        Source {
            name: "relab",
            seeds: || seeds(1, 30),
            about: "every label option on awkward headers, abundance errors, ee",
            make: relab,
        },
        Source {
            name: "sqtrim",
            seeds: || seeds(1, 20),
            about: "--fastx_filter's kept and discarded reads",
            make: sqtrim,
        },
        Source {
            name: "sqfilt",
            seeds: || seeds(1, 20),
            about: "every filter, pairs, and --fastx_mask unmasked shares",
            make: sqfilt,
        },
        Source {
            name: "rerep",
            seeds: || seeds(1, 30),
            about: "FASTA, --no-sizeout, and FASTQ through conv",
            make: rerep,
        },
        Source {
            name: "topn",
            seeds: || seeds(1, 5),
            about: "derep | topn N against --fastx_uniques --topn N",
            make: topn,
        },
        Source {
            name: "revcomp",
            seeds: || seeds(1, 20),
            about: "every nucleotide code in both cases",
            make: revcomp,
        },
        Source {
            name: "sortby",
            seeds: || seeds(1, 20),
            about: "ties of length, abundance and header; FASTQ",
            make: sortby,
        },
        Source {
            name: "getseq",
            seeds: || seeds(1, 30),
            about: "every selection mode, selected and discarded",
            make: getseq,
        },
        Source {
            name: "subseq",
            seeds: || seeds(1, 20),
            about: "many regions, FASTA and FASTQ, abundances of zero",
            make: subseq,
        },
        Source {
            name: "fqchar",
            seeds: || seeds(1, 60),
            about: "reports against vsearch's --log, several tails",
            make: fqchar,
        },
        Source {
            name: "fqstat",
            seeds: || seeds(1, 200),
            about: "reports against --log at offsets 33 and 64",
            make: fqstat,
        },
        Source {
            name: "eestat",
            seeds: || seeds(1, 200),
            about: "tables against --output at offsets 33 and 64",
            make: eestat,
        },
        Source {
            name: "eegrid",
            seeds: || seeds(1, 200),
            about: "grids against --output, several cutoffs",
            make: eegrid,
        },
        Source {
            name: "shufl",
            seeds: || seeds(1, 30),
            about: "random sizes and seeds, negative and extreme",
            make: shufl,
        },
        Source {
            name: "subsmp",
            seeds: || seeds(1, 30),
            about: "draws by count and share, by size or not",
            make: subsmp,
        },
        Source {
            name: "scrmbl",
            seeds: || seeds(1, 20),
            about: "every k from 1 to 9, random seeds",
            make: scrmbl,
        },
        Source {
            name: "lcmask",
            seeds: || seeds(1, 200),
            about: "DUST soft and hard, hard masking alone, clearing",
            make: lcmask,
        },
        Source {
            name: "syncpa",
            seeds: || seeds(1, 300),
            about: "shuffled pairs, orphans, separator sets",
            make: syncpa,
        },
        Source {
            name: "joinpa",
            seeds: || seeds(1, 300),
            about: "pairs of any lengths, five pads",
            make: joinpa,
        },
        Source {
            name: "mergpa",
            seeds: || seeds(1, 700),
            about: "pairs from random fragments, four option sets",
            make: mergpa,
        },
        Source {
            name: "redig",
            seeds: || seeds(1, 1000),
            about: "planted sites, patterns with ^ and _ anywhere",
            make: redig,
        },
    ]
}

/// vsearch's output options for a file's format: records, discarded, and
/// the width option FASTA needs.
fn fx(name: &str) -> (&'static str, &'static str, &'static str) {
    if name.ends_with(".fq") {
        ("--fastqout", "--fastqout_discarded", "")
    } else {
        ("--fastaout", "--fastaout_discarded", "--fasta_width 0")
    }
}

fn step(tool: &str, args: &[&str], stdin: Option<String>, stdout: String) -> Step {
    Step {
        tool: tool.into(),
        args: args.iter().map(|a| a.to_string()).collect(),
        stdin,
        stdout: Some(stdout),
        system: false,
    }
}

fn derep(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("derep", dir, seed, &[])?;
    let t = p(dir);
    let c = "--quiet --minseqlength 0 --notrunclabels --fasta_width 0 --sizein --sizeout";
    let mut out = Vec::new();
    for (name, f, v) in [
        ("default", "q.fq", ""),
        ("strand", "q.fq", "--strand both"),
        ("ascii64", "q64.fq", "--fastq_ascii 64"),
        (
            "ascii64-out33",
            "q64.fq",
            "--fastq_ascii 64 --fastq_asciiout 33",
        ),
        ("ascii33-out64", "q.fq", "--fastq_asciiout 64"),
        ("qout-range", "q.fq", "--fastq_qminout 5 --fastq_qmaxout 30"),
        (
            "qout-range-strand",
            "q64.fq",
            "--fastq_ascii 64 --fastq_qminout 3 --fastq_qmaxout 35 --strand both",
        ),
        ("best", "q.fq", "--fastq_qout_max"),
        (
            "best-strand",
            "q.fq",
            "--fastq_qout_max --strand both --fastq_qmaxout 39",
        ),
    ] {
        out.push(Case::new(
            name,
            argv(&format!(
                "--fastx_uniques {t}/{f} {c} --fastqout {{O}}/out {v}"
            )),
        ));
    }
    for (name, f, v) in [
        ("uc", "a.fa", ""),
        ("uc-strand", "a.fa", "--strand both"),
        ("uc-fq", "q.fq", ""),
    ] {
        let o = fx(f).0;
        out.push(Case::new(
            name,
            argv(&format!(
                "--fastx_uniques {t}/{f} {c} {o} {{O}}/out --uc {{O}}/uc --tabbedout {{O}}/tab {v}"
            )),
        ));
    }
    out.push(Case::new(
        "unsized",
        argv(&format!(
            "--fastx_uniques {t}/a.fa --quiet --minseqlength 0 --notrunclabels --fasta_width 0 \
             --strand both --fastaout {{O}}/out --uc {{O}}/uc --tabbedout {{O}}/tab"
        )),
    ));
    out.push(Case::new(
        "derep_id",
        argv(&format!(
            "--derep_id {t}/a.fa {c} --output {{O}}/out --uc {{O}}/uc"
        )),
    ));
    out.push(Case::new(
        "derep_id-strand",
        argv(&format!(
            "--derep_id {t}/a.fa {c} --strand both --output {{O}}/out"
        )),
    ));
    out.push(Case::new(
        "prefix",
        argv(&format!(
            "--derep_prefix {t}/p.fa {c} --output {{O}}/out --uc {{O}}/uc"
        )),
    ));
    out.push(Case::new(
        "prefix-nosizein",
        argv(&format!(
            "--derep_prefix {t}/p.fa --quiet --minseqlength 0 --notrunclabels --fasta_width 0 \
             --sizeout --output {{O}}/out --uc {{O}}/uc"
        )),
    ));
    out.push(Case::new(
        "smallmem",
        argv(&format!(
            "--derep_smallmem {t}/a.fa {c} --fastaout {{O}}/out"
        )),
    ));
    Ok(out)
}

fn conv(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("conv", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for (f, v) in [
        ("q33.fq", ""),
        ("q33.fq", "--fastq_asciiout 64 --fastq_qmaxout 62"),
        ("q33.fq", "--fastq_qminout 3 --fastq_qmaxout 40"),
        ("q33.fq", "--fastq_qmax 41"),
        ("q33.fq", "--fastq_qmin 2"),
        ("q64.fq", "--fastq_ascii 64"),
        ("q64.fq", "--fastq_ascii 64 --fastq_asciiout 64"),
        (
            "q64.fq",
            "--fastq_ascii 64 --fastq_qminout 2 --fastq_qmaxout 41",
        ),
        ("sol.fq", "--fastq_ascii 64 --fastq_solexa"),
        (
            "sol.fq",
            "--fastq_ascii 64 --fastq_solexa --fastq_asciiout 64",
        ),
        ("sol.fq", "--fastq_ascii 64 --fastq_solexa --fastq_qmin 0"),
        ("sol.fq", "--fastq_ascii 64"),
        ("a.fa", ""),
        ("a.fa", "--fastq_asciiout 64 --fastq_qmaxout 50"),
        ("a.fa", "--fastq_qmaxout 93"),
        ("bad.fq", ""),
        ("bad.fq", "--fastq_asciiout 64 --fastq_qmaxout 62"),
        ("bad.fa", ""),
    ] {
        let cmd = if f.ends_with(".fa") {
            "--fasta2fastq"
        } else {
            "--fastq_convert"
        };
        out.push(Case::new(
            format!("{f} {v}"),
            argv(&format!("{cmd} {t}/{f} --quiet --fastqout {{O}}/out {v}")),
        ));
    }
    // FASTA output, against --fastx_subsample, which writes every read
    for (f, v) in [
        ("q33.fq", ""),
        ("q33.fq", "--fasta_width 0"),
        ("q33.fq", "--fasta_width 7"),
        ("q64.fq", ""),
        ("a.fa", ""),
        ("a.fa", "--fasta_width 0"),
        ("a.fa", "--fasta_width 1"),
    ] {
        out.push(Case::new(
            format!("fasta {f} {v}"),
            argv(&format!(
                "--fastx_subsample {t}/{f} --quiet --sample_pct 100 --minseqlength 0 \
                 --notrunclabels --fastaout {{O}}/out {v}"
            )),
        ));
    }
    Ok(out)
}

fn sff(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("sff", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for v in [
        "",
        "--sff_clip",
        "--fastq_asciiout 64 --fastq_qmaxout 50 --fastq_qminout 3",
    ] {
        out.push(Case::new(
            format!("sff {v}"),
            argv(&format!(
                "--sff_convert {t}/r.sff --quiet --fastqout {{O}}/out {v}"
            )),
        ));
    }
    // FASTA: the bases --sff_convert writes, against conv --fasta
    for clip in [false, true] {
        let sff = format!("{t}/r.sff");
        let probe = Box::new(move |rdir: &Path| {
            let o = p(&rdir.join("fasta"));
            let mut args = vec!["-i", sff.as_str(), "--fasta"];
            if clip {
                args.push("--clip");
            }
            (
                vec![step("conv", &args, None, o.clone())],
                vec![("fastqout".to_string(), o)],
            )
        });
        let flag = if clip { " --sff_clip" } else { "" };
        out.push(
            Case::new(
                format!("fasta clip={clip}"),
                argv(&format!(
                    "--sff_convert {t}/r.sff --quiet --fastqout {{O}}/out{flag}"
                )),
            )
            .probe(probe)
            .comparator("fastqout", "fastq-to-fasta"),
        );
    }
    let data = std::fs::read(dir.join("r.sff"))?;
    for cut in 1..=3u64 {
        let n = Rng::new(seed * 10 + cut).randint(0, data.len() as i64) as usize;
        std::fs::write(dir.join(format!("t{cut}.sff")), &data[..n])?;
        out.push(Case::new(
            format!("cut at {n}"),
            argv(&format!(
                "--sff_convert {t}/t{cut}.sff --quiet --fastqout {{O}}/out"
            )),
        ));
    }
    Ok(out)
}

fn relab(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("relab", dir, seed, &[])?;
    let t = p(dir);
    let opts = [
        "",
        "--relabel R",
        "--relabel '' --relabel_keep --sizeout",
        "--relabel_self --relabel_keep",
        "--relabel_sha1 --lengthout",
        "--relabel_md5 --sizeout --sample S",
        "--sizeout",
        "--xsize",
        "--xee",
        "--lengthout",
        "--xlength",
        "--xsize --xee --xlength",
        "--sample S1 --label_suffix ';'",
        "--label_suffix ';k=v;' --sizeout --lengthout",
        "--label_suffix _s --xsize --sizeout",
        "--relabel_keep --xsize",
    ];
    let mut out = Vec::new();
    let mut c = |f: &str, v: &str| {
        let (o, cmd, w) = if f.ends_with(".fa") {
            ("--fastaout", "--fastx_filter", "--fasta_width 0")
        } else if f == "e.fq" {
            ("--fastqout", "--fastq_convert", "")
        } else {
            ("--fastqout", "--fastx_filter", "")
        };
        out.push(Case::new(
            format!("{f} {v}"),
            argv(&format!("{cmd} {t}/{f} --quiet {o} {{O}}/out {w} {v}")),
        ));
    };
    for f in ["q33.fq", "a.fa", "e.fq", "bad.fa"] {
        for v in opts {
            c(f, v);
        }
    }
    for (f, v) in [
        ("q33.fq", "--eeout"),
        ("q33.fq", "--eeout --sizeout --lengthout --sample S"),
        ("q33.fq", "--eeout --relabel R --label_suffix ';'"),
        ("q33.fq", "--eeout --xee"),
        ("q64.fq", "--eeout --fastq_ascii 64"),
        ("q64.fq", "--eeout --fastq_ascii 64 --xsize --lengthout"),
        ("q33.fq", "--eeout --fastq_ascii 64"),
    ] {
        c(f, v);
    }
    Ok(out)
}

fn kept_discarded(t: &str, f: &str, cmd: &str, v: &str) -> Case {
    let (o, d, w) = fx(f);
    Case::new(
        format!("{f} {v}"),
        argv(&format!(
            "{cmd} {t}/{f} --quiet --notrunclabels {o} {{O}}/k {d} {{O}}/x {w} {v}"
        )),
    )
}

fn sqtrim(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("sqtrim", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    let mut c = |f: &str, v: &str| {
        out.push(kept_discarded(&t, f, "--fastx_filter", v).success_outputs_only())
    };
    for f in ["q33.fq", "a.fa"] {
        for v in [
            "",
            "--fastq_stripleft 3",
            "--fastq_stripright 7",
            "--fastq_stripleft 20 --fastq_stripright 30",
            "--fastq_trunclen_keep 25",
            "--fastq_stripleft 2 --fastq_trunclen_keep 10",
        ] {
            c(f, v);
        }
    }
    for q in [0, 2, 15, 25, 30] {
        c("q33.fq", &format!("--fastq_truncqual {q}"));
    }
    for e in ["0", "0.01", "0.1", "0.5", "1", "2.5"] {
        c("q33.fq", &format!("--fastq_truncee {e}"));
    }
    for e in ["0", "0.001", "0.005", "0.02", "0.1"] {
        c("q33.fq", &format!("--fastq_truncee_rate {e}"));
    }
    c(
        "q33.fq",
        "--fastq_stripleft 1 --fastq_stripright 2 --fastq_trunclen_keep 30 --fastq_truncqual 12 \
         --fastq_truncee 0.8 --fastq_truncee_rate 0.03",
    );
    c(
        "q64.fq",
        "--fastq_ascii 64 --fastq_truncqual 15 --fastq_truncee 1",
    );
    c("q64.fq", "--fastq_truncqual 15");
    c("hi.fq", "--fastq_qmax 41");
    c("hi.fq", "--fastq_qmax 41 --fastq_truncqual 20");
    c("hi.fq", "--fastq_qmax 41 --fastq_trunclen_keep 5");
    c("hi.fq", "--fastq_qmax 60 --fastq_truncee 0.05");
    c("hi.fq", "--fastq_qmin 5 --fastq_stripleft 4");
    c("bad.fq", "");
    c("bad.fq", "--fastq_stripleft 2 --fastq_truncqual 15");
    Ok(out)
}

fn sqfilt(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("sqfilt", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    // rsearch exits 2 on invalid arguments, as clap does, where vsearch exits 1
    let kw = |c: Case| c.success_outputs_only().exit_map(1, 2, 1);
    for f in ["q33.fq", "a.fa"] {
        for v in [
            "",
            "--fastq_minlen 20",
            "--fastq_maxlen 25",
            "--fastq_maxns 0",
            "--fastq_maxns 2",
            "--minsize 3 --maxsize 7",
            "--fastq_minlen 10 --fastq_maxlen 30 --fastq_maxns 1 --minsize 2",
        ] {
            out.push(kw(kept_discarded(&t, f, "--fastx_filter", v)));
        }
    }
    let mut c = |f: &str, v: &str| out.push(kw(kept_discarded(&t, f, "--fastx_filter", v)));
    for e in ["0.01", "0.1", "0.5", "1", "3"] {
        c("q33.fq", &format!("--fastq_maxee {e}"));
    }
    for e in ["0", "0.001", "0.005", "0.02"] {
        c("q33.fq", &format!("--fastq_maxee_rate {e}"));
    }
    for q in [0, 5, 15, 25] {
        c("q33.fq", &format!("--fastq_minqual {q}"));
    }
    c(
        "q64.fq",
        "--fastq_ascii 64 --fastq_maxee 0.5 --fastq_minqual 10",
    );
    c("q64.fq", "--fastq_maxee 0.5");
    c("hi.fq", "--fastq_qmax 41");
    c("hi.fq", "--fastq_qmax 60 --fastq_maxee 0.01");
    c("hi.fq", "--fastq_qmin 3");
    for v in [
        "",
        "--fastq_maxee 0.5",
        "--fastq_minlen 15 --fastq_maxns 1",
        "--fastq_minqual 10 --minsize 3",
    ] {
        out.push(kw(Case::new(
            format!("pairs {v}"),
            argv(&format!(
                "--fastx_filter {t}/f.fq --reverse {t}/r.fq --quiet --notrunclabels --fastqout {{O}}/k \
                 --fastqout_discarded {{O}}/x --fastqout_rev {{O}}/k2 --fastqout_discarded_rev {{O}}/x2 {v}"
            )),
        )));
    }
    let mut m = |f: &str, v: &str| {
        let o = if f.ends_with(".fq") {
            "--fastqout"
        } else {
            "--fasta_width 0 --fastaout"
        };
        out.push(kw(Case::new(
            format!("mask {f} {v}"),
            argv(&format!(
                "--fastx_mask {t}/{f} --quiet --notrunclabels --qmask soft {o} {{O}}/k {v}"
            )),
        )));
    };
    for x in ["0", "33.3", "50", "80", "100"] {
        m("a.fa", &format!("--min_unmasked_pct {x}"));
        m("a.fa", &format!("--max_unmasked_pct {x}"));
        m(
            "q33.fq",
            &format!("--min_unmasked_pct {x} --max_unmasked_pct 95"),
        );
        m("h.fa", &format!("--hardmask --min_unmasked_pct {x}"));
        m("h.fa", &format!("--hardmask --max_unmasked_pct {x}"));
    }
    Ok(out)
}

fn rerep(dir: &Path, seed: u64, env: &Env) -> io::Result<Vec<Case>> {
    generate("rerep", dir, seed, &[])?;
    let t = p(dir);
    let c = "--quiet --notrunclabels --fasta_width 0 --output {O}/out";
    let w = ["missing-abundance"];
    let mut out = vec![
        Case::new(
            "sizeout",
            argv(&format!("--rereplicate {t}/a.fa {c} --sizeout")),
        )
        .warn(&w),
        Case::new("no-sizeout", argv(&format!("--rereplicate {t}/a.fa {c}"))).warn(&w),
    ];
    // FASTQ, which vsearch does not read: rerep | conv --fasta against
    // --rereplicate on the same reads in FASTA
    let _ = std::process::Command::new(env.bindir.join("conv"))
        .args([
            "-i",
            &format!("{t}/q.fq"),
            "--fasta",
            "--fasta-width",
            "0",
            "-o",
            &format!("{t}/q.fa"),
        ])
        .stderr(std::process::Stdio::null())
        .status();
    let q = format!("{t}/q.fq");
    let probe = Box::new(move |rdir: &Path| {
        let (tmp, o) = (p(&rdir.join("tmp1")), p(&rdir.join("out")));
        (
            vec![
                step("rerep", &["-i", q.as_str()], None, tmp.clone()),
                step(
                    "conv",
                    &["--fasta", "--fasta-width", "0"],
                    Some(tmp),
                    o.clone(),
                ),
            ],
            vec![("output".to_string(), o)],
        )
    });
    out.push(
        Case::new(
            "fastq",
            argv(&format!("--rereplicate {t}/q.fa {c} --sizeout")),
        )
        .warn(&w)
        .probe(probe),
    );
    Ok(out)
}

fn topn(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("derep", dir, seed, &[])?;
    let t = p(dir);
    let c = "--quiet --sizein --sizeout --notrunclabels --minseqlength 0 --fasta_width 0";
    let mut out = Vec::new();
    for f in ["a.fa", "q.fq"] {
        let o = fx(f).0;
        for n in [1, 2, 3, 5, 8, 13, 50, 1000] {
            for s in ["", "--strand both"] {
                out.push(Case::new(
                    format!("{f} topn {n} {s}"),
                    argv(&format!(
                        "--fastx_uniques {t}/{f} {c} {o} {{O}}/out --topn {n} {s}"
                    )),
                ));
            }
        }
    }
    Ok(out)
}

fn revcomp(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("revcomp", dir, seed, &[])?;
    let t = p(dir);
    Ok(["a.fa", "q.fq"]
        .iter()
        .map(|f| {
            Case::new(
                *f,
                argv(&format!(
                    "--fastx_revcomp {t}/{f} --quiet --notrunclabels --fasta_width 0 {} {{O}}/out",
                    fx(f).0
                )),
            )
        })
        .collect())
}

fn sortby(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("sortby", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for k in ["length", "size"] {
        for f in ["a.fa", "q.fq"] {
            out.push(Case::new(
                format!("{k} {f}"),
                argv(&format!(
                    "--sortby{k} {t}/{f} --quiet --notrunclabels --minseqlength 0 --fasta_width 0 --output {{O}}/out"
                )),
            ));
        }
    }
    Ok(out)
}

fn getseq(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("getseq", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    let mut c = |cmd: &str, f: &str, v: Vec<String>| {
        let (o, nm, w) = if f.ends_with(".fq") {
            ("--fastqout", "--notmatchedfq", "")
        } else {
            ("--fastaout", "--notmatched", "--fasta_width 0")
        };
        let mut a = argv(&format!(
            "{cmd} {t}/{f} --quiet --notrunclabels {w} {o} {{O}}/k {nm} {{O}}/x"
        ));
        let name = format!("{f} {cmd} {}", v.join(" "));
        a.extend(v);
        out.push(Case::new(name, a));
    };
    let s = |w: &[&str]| w.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    for f in ["a.fa", "q.fq"] {
        for l in ["seq1", "SEQ1", "s1", "x", "ok", "a b", "seq1;sample=s1"] {
            c("--fastx_getseq", f, s(&["--label", l]));
            c(
                "--fastx_getseq",
                f,
                s(&["--label", l, "--label_substr_match"]),
            );
            c("--fastx_getseqs", f, s(&["--label_word", l]));
            c(
                "--fastx_getseqs",
                f,
                s(&["--label_word", l, "--label_field", "sample"]),
            );
            c(
                "--fastx_getseqs",
                f,
                s(&["--label_word", l, "--label_field", "tag"]),
            );
        }
        let labels = format!("{t}/labels.txt");
        let words = format!("{t}/words.txt");
        c("--fastx_getseqs", f, s(&["--labels", &labels]));
        c(
            "--fastx_getseqs",
            f,
            s(&["--labels", &labels, "--label_substr_match"]),
        );
        c("--fastx_getseqs", f, s(&["--label_words", &words]));
        c(
            "--fastx_getseqs",
            f,
            s(&["--label_words", &words, "--label_field", "sample"]),
        );
    }
    for l in ["seq1", "s1"] {
        c("--fastx_getseq", "bad.fa", s(&["--label", l]));
        c("--fastx_getseqs", "bad.fa", s(&["--label_word", l]));
    }
    Ok(out)
}

fn subseq(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("subseq", dir, seed, &[])?;
    let t = p(dir);
    let mut r = Rng::new(seed);
    let regions: Vec<(i64, i64)> = (0..8)
        .map(|_| {
            let s = r.randint(1, 45);
            (s, r.randint(s, 60))
        })
        .collect();
    let mut out = Vec::new();
    for f in ["a.fa", "q.fq", "bad.fa", "bad.fq"] {
        let (o, w) = if f.ends_with(".fq") {
            ("--fastqout", "")
        } else {
            ("--fastaout", "--fasta_width 0")
        };
        let mut base = argv(&format!(
            "--fastx_getsubseq {t}/{f} --quiet --notrunclabels --label_substr_match {w} {o} {{O}}/out"
        ));
        base.extend(["--label".to_string(), String::new()]);
        out.push(Case::new(f, base.clone()));
        for (s, e) in &regions {
            for v in [
                vec![
                    "--subseq_start".into(),
                    s.to_string(),
                    "--subseq_end".into(),
                    e.to_string(),
                ],
                vec!["--subseq_start".into(), s.to_string()],
                vec!["--subseq_end".into(), e.to_string()],
            ] {
                let mut a = base.clone();
                let name = format!("{f} {}", v.join(" "));
                a.extend(v);
                out.push(Case::new(name, a));
            }
        }
    }
    Ok(out)
}

fn fqchar(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("fqchar", dir, seed, &[])?;
    let t = p(dir);
    Ok([1, 2, 4, 7]
        .iter()
        .map(|n| {
            Case::new(
                format!("tail {n}"),
                argv(&format!(
                    "--fastq_chars {t}/q.fq --quiet --fastq_tail {n} --log {{O}}/log"
                )),
            )
        })
        .collect())
}

const STAT_CASES: [(&str, &str); 4] = [
    ("q.fq", ""),
    ("q.fq", "--fastq_qmax 41"),
    ("s.fq", "--fastq_ascii 64 --fastq_qmin -5"),
    ("s.fq", "--fastq_ascii 64"),
];

fn fqstat(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("fqstat", dir, seed, &[])?;
    let t = p(dir);
    Ok(STAT_CASES
        .iter()
        .map(|(f, v)| {
            Case::new(
                format!("{f} {v}"),
                argv(&format!(
                    "--fastq_stats {t}/{f} --quiet --log {{O}}/log {v}"
                )),
            )
            .success_outputs_only()
        })
        .collect())
}

fn report_cases(t: &str, cmd: &str, cases: &[(&str, &str)]) -> Vec<Case> {
    cases
        .iter()
        .map(|(f, v)| {
            Case::new(
                format!("{f} {v}"),
                argv(&format!("{cmd} {t}/{f} --quiet --output {{O}}/out {v}")),
            )
            .success_outputs_only()
        })
        .collect()
}

fn eestat(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("fqstat", dir, seed, &[])?;
    Ok(report_cases(&p(dir), "--fastq_eestats", &STAT_CASES))
}

fn eegrid(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("fqstat", dir, seed, &[])?;
    Ok(report_cases(
        &p(dir),
        "--fastq_eestats2",
        &[
            ("q.fq", ""),
            ("q.fq", "--length_cutoffs 1,*,1 --ee_cutoffs 0.1,1,0.5,3"),
            ("q.fq", "--length_cutoffs 3,100,7 --ee_cutoffs 0.25"),
            ("q.fq", "--length_cutoffs 200,*,25 --fastq_qmax 41"),
            (
                "s.fq",
                "--length_cutoffs 2,*,3 --ee_cutoffs 1,2,5 --fastq_ascii 64 --fastq_qmin -5",
            ),
            ("s.fq", "--fastq_ascii 64"),
        ],
    ))
}

fn shufl(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("shufl", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for s in read(&dir.join("seeds")).split_whitespace() {
        for f in ["a.fa", "q.fq"] {
            out.push(Case::new(
                format!("{f} seed {s}"),
                argv(&format!(
                    "--shuffle {t}/{f} --quiet --notrunclabels --minseqlength 0 --fasta_width 0 \
                     --randseed {s} --output {{O}}/out"
                )),
            ));
        }
    }
    Ok(out)
}

fn subsmp(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("subsmp", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for line in read(&dir.join("draws")).lines().filter(|l| !l.is_empty()) {
        let Some((s, size)) = line.split_once('|') else {
            continue;
        };
        let vsize = size
            .replace("--count", "--sample_size")
            .replace("--pct", "--sample_pct");
        for f in ["a.fa", "q.fq"] {
            let (o, d, w) = fx(f);
            for by in ["", "--sizein --sizeout"] {
                out.push(Case::new(
                    format!("{f} seed {s} {size} {by}"),
                    argv(&format!(
                        "--fastx_subsample {t}/{f} --quiet --notrunclabels --minseqlength 0 {w} {o} \
                         {{O}}/k {d} {{O}}/x --randseed {s} {vsize} {by}"
                    )),
                ));
            }
        }
    }
    Ok(out)
}

fn scrmbl(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("scrmbl", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for s in read(&dir.join("seeds")).split_whitespace() {
        for k in 1..=9 {
            for f in ["a.fa", "q.fq"] {
                let (o, _, w) = fx(f);
                out.push(Case::new(
                    format!("{f} seed {s} k {k}"),
                    argv(&format!(
                        "--scramble {t}/{f} --quiet --notrunclabels {w} {o} {{O}}/out --randseed {s} \
                         --scramble_kmer {k}"
                    )),
                ));
            }
        }
    }
    Ok(out)
}

fn lcmask(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("lcmask", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for f in ["a.fa", "q.fq"] {
        let (o, _, w) = fx(f);
        for v in [
            "--qmask dust",
            "--qmask dust --hardmask",
            "--qmask soft --hardmask",
        ] {
            out.push(Case::new(
                format!("{f} {v}"),
                argv(&format!(
                    "--fastx_mask {t}/{f} --quiet --notrunclabels --minseqlength 0 {w} {o} {{O}}/out {v}"
                )),
            ));
        }
        // clearing: vsearch writes the input as it is, upper-cased here on
        // the sequence lines
        let input = format!("{t}/{f}");
        let key = o.trim_start_matches('-').to_string();
        let probe_key = key.clone();
        let probe = Box::new(move |rdir: &Path| {
            let out = p(&rdir.join("out"));
            (
                vec![step(
                    "lcmask",
                    &["-i", input.as_str(), "--clear"],
                    None,
                    out.clone(),
                )],
                vec![(probe_key.clone(), out)],
            )
        });
        out.push(
            Case::new(
                format!("{f} --clear"),
                argv(&format!(
                    "--fastx_mask {t}/{f} --quiet --notrunclabels --minseqlength 0 --fasta_width 0 \
                     --qmask none {o} {{O}}/out"
                )),
            )
            .probe(probe)
            .comparator(&key, "upper-seq"),
        );
    }
    Ok(out)
}

fn syncpa(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("syncpa", dir, seed, &[])?;
    let t = p(dir);
    let mut out = Vec::new();
    for sep in ["", "_/", "._/"] {
        for (f, r, x) in [("f.fa", "r.fa", "fasta"), ("f.fq", "r.fq", "fastq")] {
            let mut v = argv(&format!(
                "--fastx_syncpairs {t}/{f} --reverse {t}/{r} --quiet --fasta_width 0 --{x}out {{O}}/1 \
                 --{x}out_rev {{O}}/2 --{x}out_orphans {{O}}/3 --{x}out_orphans_rev {{O}}/4"
            ));
            if !sep.is_empty() {
                v.extend(["--read_separators".to_string(), sep.to_string()]);
            }
            out.push(Case::new(format!("{f} sep[{sep}]"), v).errors("exact"));
        }
    }
    Ok(out)
}

fn joinpa(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("joinpa", dir, seed, &[])?;
    let t = p(dir);
    let base = argv(&format!(
        "--fastq_join {t}/f.fq --reverse {t}/r.fq --quiet --fastqout {{O}}/out"
    ));
    let pads: [&[&str]; 5] = [
        &[],
        &["--join_padgap", "", "--join_padgapq", ""],
        &["--join_padgap", "ACG"],
        &["--join_padgap", "x y", "--join_padgapq", "!~#"],
        &["--join_padgap", "NN", "--fastq_ascii", "64"],
    ];
    Ok(pads
        .iter()
        .enumerate()
        .map(|(k, v)| {
            let mut a = base.clone();
            a.extend(v.iter().map(|x| x.to_string()));
            Case::new(format!("pad {k}"), a).errors("joinpa")
        })
        .collect())
}

fn mergpa(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("mergpa", dir, seed, &[])?;
    let t = p(dir);
    let base = format!(
        "--fastq_mergepairs {t}/f.fq --reverse {t}/r.fq --quiet --fastqout {{O}}/1 \
         --fastqout_notmerged_fwd {{O}}/2 --fastqout_notmerged_rev {{O}}/3 --eetabbedout {{O}}/4 --log {{O}}/log"
    );
    Ok([
        "",
        "--fastq_minovlen 5 --fastq_maxdiffs 3 --fastq_allowmergestagger",
        "--fastq_minovlen 7 --fastq_maxdiffpct 5 --fastq_qmaxout 60 --fastq_qminout 3",
        "--fastq_minovlen 20 --fastq_maxdiffs 0 --fastq_qmin 2 --fastq_qmax 44",
    ]
    .iter()
    .enumerate()
    .map(|(k, v)| {
        Case::new(format!("opts {k}"), argv(&format!("{base} {v}")))
            .errors("mergpa")
            .success_outputs_only()
    })
    .collect())
}

fn redig(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("redig", dir, seed, &[])?;
    let t = p(dir);
    let pattern = read(&dir.join("pattern"));
    let mut a = argv(&format!(
        "--cut {t}/in.fa --quiet --fasta_width 0 --fastaout {{O}}/1 --fastaout_rev {{O}}/2 \
         --fastaout_discarded {{O}}/3 --fastaout_discarded_rev {{O}}/4"
    ));
    a.extend(["--cut_pattern".to_string(), pattern.clone()]);
    Ok(vec![
        Case::new(format!("pattern {pattern}"), a).errors("abundance"),
    ])
}
