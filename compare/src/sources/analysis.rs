//! Cases for the analysis tools, from the former zsh runners with the
//! same generators and seeds.
//!
//! Where a runner applied an rsearch tool to vsearch's own hits (hitfmt,
//! otutab, lcatax, getseq --hits, clspl), the case is vsearch's command
//! whole. The translated pipeline runs the search too, which checks the
//! search and the formatting together.

use std::io;
use std::path::Path;
use std::process::{Command, Stdio};

use super::{Env, Source, argv, generate, p, params, read, seeds, strs};
use crate::cases::Case;
use crate::execute::Run;

pub fn sources() -> Vec<Source> {
    vec![
        Source {
            name: "mkudb",
            seeds: || seeds(1, 600),
            about: "UDB files through sqfilt and lcmask, six masking choices, word lengths 3 to 13",
            make: mkudb,
        },
        Source {
            name: "mkudb-long",
            seeds: || vec![13, 14, 15],
            about: "the long words: files of 4^k words (256 MB at 13, 4 GB at 15)",
            make: mkudb_long,
        },
        Source {
            name: "udbinfo",
            seeds: || seeds(1, 500),
            about: "headers valid or with one field wrong",
            make: udbinfo,
        },
        Source {
            name: "udbstat",
            seeds: || seeds(1, 700),
            about: "mkudb's files, half of them spoiled in one part; vsearch aborts skipped",
            make: udbstat,
        },
        Source {
            name: "conv-udb",
            seeds: || seeds(1, 700),
            about: "udbstat's spoiled files against --udb2fasta, unwrapped and at 80",
            make: conv_udb,
        },
        Source {
            name: "dbornt",
            seeds: || seeds(1, 600),
            about: "queries cut from a database, six database masking choices, three query ones, FASTA or UDB",
            make: dbornt,
        },
        Source {
            name: "sintax",
            seeds: || seeds(1, 600),
            about: "taxonomies on random trees, strands, cutoffs, random ties, seeds, masking, FASTA or UDB",
            make: sintax,
        },
        Source {
            name: "exsrch",
            seeds: || seeds(1, 600),
            about: "few sequences copied in either case, strands, --maxhits, --self, abundance filters",
            make: exsrch,
        },
        Source {
            name: "hitfmt",
            seeds: || seeds(1, 600),
            about: "--usearch_global's hit outputs: every identity definition, scores, row lengths, DUST masking",
            make: hitfmt,
        },
        Source {
            name: "otutab",
            seeds: || seeds(1, 600),
            about: "the tables of exact and global searches",
            make: otutab,
        },
        Source {
            name: "lcatax",
            seeds: || seeds(1, 600),
            about: "taxonomies on random tree paths, many hits, cutoffs",
            make: lcatax,
        },
        Source {
            name: "getseq-hits",
            seeds: || seeds(1, 600),
            about: "--usearch_global's matched and not matched queries and targets, from hitfmt's searches",
            make: getseq_hits,
        },
        Source {
            name: "allpr",
            seeds: || seeds(1, 1300),
            about: "every criterion, identity definition, score and gap penalty",
            make: allpr,
        },
        Source {
            name: "dbsrch",
            seeds: || seeds(1, 1000),
            about: "related sequences, UDB databases, accepts and rejects, words, strands, criteria, scores",
            make: dbsrch,
        },
        Source {
            name: "clust",
            seeds: || seeds(1, 1000),
            about: "input order or sorted by length or size; limits, words, strands, abundances, criteria, scores",
            make: clust,
        },
        Source {
            name: "unoise",
            seeds: || seeds(1, 1000),
            about: "abundant templates and their variants; alpha, identity, limits, words, strands, abundances",
            make: unoise,
        },
        Source {
            name: "clspl",
            seeds: || seeds(1, 600),
            about: "--clusters, from clust's random inputs",
            make: clspl,
        },
        Source {
            name: "clmsa",
            seeds: || seeds(1, 600),
            about: "--msaout, --consout and --profile, from clust's random inputs",
            make: clmsa,
        },
        Source {
            name: "uchime",
            seeds: || seeds(1, 1000),
            about: "the four UCHIME commands, scores, thresholds, skew, alignment width; all five outputs",
            make: uchime,
        },
        Source {
            name: "mchime",
            seeds: || seeds(1, 1000),
            about: "chimeras of two to five long parents; parts, parents, segments, differences, skew, width",
            make: mchime,
        },
    ]
}

const VFLAGS: &[&str] = &[
    "acceptall",
    "leftjust",
    "rightjust",
    "self",
    "selfid",
    "top_hits_only",
    "n_mismatch",
    "sizeorder",
];

/// A generator's parameter words (vsearch option names without dashes,
/// values after them) as vsearch options.
fn vopts(words: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < words.len() {
        let n = &words[i];
        out.push(format!("--{n}"));
        if !VFLAGS.contains(&n.as_str()) {
            if let Some(v) = words.get(i + 1) {
                out.push(v.clone());
            }
            i += 1;
        }
        i += 1;
    }
    out
}

/// Runs an rsearch tool (to make a UDB database both programs read).
fn sh(env: &Env, tool: &str, args: &[&str]) -> bool {
    Command::new(env.bindir.join(tool))
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

/// Runs rsearch tools one into the next, from src into dst.
fn pipe(env: &Env, cmds: &[(&str, Vec<&str>)], src: &str, dst: &str) -> bool {
    let mut cur = src.to_string();
    for (k, (tool, args)) in cmds.iter().enumerate() {
        let out = if k + 1 == cmds.len() {
            dst.to_string()
        } else {
            format!("{dst}.{k}")
        };
        let Ok(f) = std::fs::File::create(&out) else {
            return false;
        };
        let ok = Command::new(env.bindir.join(tool))
            .arg("-i")
            .arg(&cur)
            .args(args)
            .stdout(f)
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|s| s.success());
        if !ok {
            return false;
        }
        cur = out;
    }
    true
}

fn nonempty(path: &str) -> bool {
    std::fs::metadata(path).is_ok_and(|m| m.len() > 0)
}

fn mask_v(mask: &str) -> Vec<String> {
    match mask {
        "dusthard" => strs(&["--dbmask", "dust", "--hardmask"]),
        "softhard" => strs(&["--dbmask", "soft", "--hardmask"]),
        "nonehard" => strs(&["--dbmask", "none", "--hardmask"]),
        m => strs(&["--dbmask", m]),
    }
}

fn remove(dir: &Path, names: &[&str]) {
    for n in names {
        let _ = std::fs::remove_file(dir.join(n));
    }
}

fn mkudb_case(dir: &Path) -> Case {
    let t = p(dir);
    let w = params(dir);
    let (k, minlen, maxlen, mask, fmt) = (&w[0], &w[1], &w[2], &w[3], &w[4]);
    let mut a = argv(&format!(
        "--makeudb_usearch {t}/in.{fmt} --output {{O}}/udb --quiet --notrunclabels --wordlength {k} \
         --minseqlength {minlen} --maxseqlength {maxlen}"
    ));
    a.extend(mask_v(mask));
    Case::new(format!("k={k} mask={mask} fmt={fmt}"), a)
        .errors("mkudb")
        .warn(&["nowords"])
        .success_only(&["output", "nowords"])
        .cover("k", k.as_str())
        .cover("mask", mask.as_str())
        .cover("fmt", fmt.as_str())
}

fn mkudb(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    remove(dir, &["in.fa", "in.fq"]);
    generate("mkudb", dir, seed, &[])?;
    Ok(vec![mkudb_case(dir)])
}

fn mkudb_long(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    remove(dir, &["in.fa", "in.fq"]);
    let k = seed.to_string();
    generate("mkudb", dir, seed, &[&k])?;
    Ok(vec![mkudb_case(dir)])
}

fn udbinfo(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("udbinfo", dir, seed, &[])?;
    let t = p(dir);
    Ok(vec![
        Case::new(
            "udbinfo",
            argv(&format!("--udbinfo {t}/in.udb --quiet --log {{O}}/log")),
        )
        .errors("invalid"),
    ])
}

fn udb_skip(udb: String, magic: bool) -> crate::cases::Skip {
    Box::new(move |v: Option<&Run>, _: Option<&Run>| {
        let v = v?;
        if v.exit == 134 {
            return Some("vsearch aborted".into());
        }
        if magic {
            let head = std::fs::read(&udb).unwrap_or_default();
            if !head.starts_with(b"FBDU") && v.stderr.contains("Invalid UDB file") {
                return Some("not a UDB file to conv".into());
            }
        }
        None
    })
}

fn udbstat(dir: &Path, seed: u64, env: &Env) -> io::Result<Vec<Case>> {
    let mkudb = p(&env.bindir.join("mkudb"));
    generate("udbstat", dir, seed, &[&mkudb])?;
    let t = p(dir);
    Ok(vec![
        Case::new(
            "udbstats",
            argv(&format!("--udbstats {t}/in.udb --quiet --log {{O}}/log")),
        )
        .errors("udb")
        .success_only(&["log"])
        .skip(udb_skip(format!("{t}/in.udb"), false))
        .cover("what", read(&dir.join("what")).trim()),
    ])
}

fn conv_udb(dir: &Path, seed: u64, env: &Env) -> io::Result<Vec<Case>> {
    let mkudb = p(&env.bindir.join("mkudb"));
    generate("udbstat", dir, seed, &[&mkudb])?;
    let t = p(dir);
    Ok([0, 80]
        .iter()
        .map(|w| {
            Case::new(
                format!("width {w}"),
                argv(&format!(
                    "--udb2fasta {t}/in.udb --output {{O}}/fa --quiet --fasta_width {w}"
                )),
            )
            .errors("udb")
            .skip(udb_skip(format!("{t}/in.udb"), true))
        })
        .collect())
}

fn dbornt(dir: &Path, seed: u64, env: &Env) -> io::Result<Vec<Case>> {
    remove(dir, &["q.fa", "q.fq"]);
    generate("dbornt", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let (k, dbmask, qmask, dbkind, fmt) = (&w[0], &w[1], &w[2], &w[3], &w[4]);
    let (mut db, mut vmask) = (format!("{t}/db.fa"), mask_v(dbmask));
    if dbkind == "udb" {
        let masking: Vec<(&str, Vec<&str>)> = match dbmask.as_str() {
            "dust" => vec![("lcmask", vec![])],
            "dusthard" => vec![("lcmask", vec!["--hard"])],
            "soft" => vec![],
            "softhard" => vec![("lcmask", vec!["--hard", "--no-dust"])],
            _ => vec![("lcmask", vec!["--clear"])],
        };
        let mut cmds = vec![("sqfilt", vec!["--max-len", "50000"])];
        cmds.extend(masking);
        let rdb = format!("{t}/rdb.fa");
        let udb = format!("{t}/db.udb");
        if pipe(env, &cmds, &db, &rdb)
            && nonempty(&rdb)
            && sh(env, "mkudb", &["-i", &rdb, "--word-len", k, "-o", &udb])
        {
            db = udb;
            vmask = vec![];
        }
    }
    let out = if fmt == "fq" {
        "--fastqout"
    } else {
        "--fastaout"
    };
    let records = if qmask == "none" {
        String::new()
    } else {
        format!("{out} {{O}}/out --notmatched {{O}}/not")
    };
    let mut a = argv(&format!(
        "--orient {t}/q.{fmt} --db {db} {records} --tabbedout {{O}}/tab --quiet --notrunclabels \
         --fasta_width 0 --wordlength {k} --qmask {qmask} --log {{O}}/log"
    ));
    a.extend(vmask);
    Ok(vec![
        Case::new(
            format!("k={k} dbmask={dbmask} qmask={qmask} db={dbkind} fmt={fmt}"),
            a,
        )
        .errors("abundance")
        .warn(&["nowords", "wordlength-adjusted"])
        .success_only(&["log", "nowords", "wordlength-adjusted"])
        .cover("dbmask", dbmask.as_str())
        .cover("qmask", qmask.as_str())
        .cover("db", dbkind.as_str())
        .count("fwd", "tabbedout", "col 2 +")
        .count("rev", "tabbedout", "col 2 -")
        .count("und", "tabbedout", "col 2 ?"),
    ])
}

fn sintax(dir: &Path, seed: u64, env: &Env) -> io::Result<Vec<Case>> {
    remove(dir, &["q.fa", "q.fq"]);
    generate("sintax", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let (k, strand, cutoff, rand, rseed, dbmask, qmode, dbkind, fmt) = (
        &w[0], &w[1], &w[2], &w[3], &w[4], &w[5], &w[6], &w[7], &w[8],
    );
    let mut db = format!("{t}/db.fa");
    let mut q = format!("{t}/q.{fmt}");
    if dbkind == "udb" {
        let rdb = format!("{t}/rdb.fa");
        let udb = format!("{t}/db.udb");
        let mut cmds = vec![("sqfilt", vec!["--min-len", "32", "--max-len", "50000"])];
        if dbmask == "none" {
            cmds.push(("lcmask", vec!["--clear"]));
        }
        if pipe(env, &cmds, &db, &rdb)
            && nonempty(&rdb)
            && sh(env, "mkudb", &["-i", &rdb, "--word-len", k, "-o", &udb])
        {
            db = udb;
        }
    }
    if qmode == "mask" {
        // vsearch never masks queries: their lower case turned into N
        let text = read(Path::new(&q));
        let lines: Vec<String> = text
            .split('\n')
            .enumerate()
            .map(|(i, l)| {
                let seq = if fmt == "fq" {
                    i % 4 == 1
                } else {
                    !l.starts_with('>')
                };
                if seq {
                    l.chars()
                        .map(|c| if c.is_lowercase() { 'N' } else { c })
                        .collect()
                } else {
                    l.into()
                }
            })
            .collect();
        q = format!("{t}/v.q");
        std::fs::write(&q, lines.join("\n"))?;
    }
    let mut a = argv(&format!(
        "--sintax {q} --db {db} --tabbedout {{O}}/tab --quiet --notrunclabels --threads 1 --strand {strand} \
         --randseed {rseed} --wordlength {k} --dbmask {dbmask} --log {{O}}/log"
    ));
    if cutoff != "0" {
        a.extend(strs(&["--sintax_cutoff", cutoff]));
    }
    if rand == "1" {
        a.push("--sintax_random".into());
    }
    Ok(vec![
        Case::new(
            format!(
                "k={k} strand={strand} cutoff={cutoff} random={rand} dbmask={dbmask} q={qmode} db={dbkind} fmt={fmt}"
            ),
            a,
        )
        .errors("abundance")
        .success_only(&["log"])
        .cover("strand", strand.as_str())
        .cover("db", dbkind.as_str())
        .cover("q", qmode.as_str())
        .count("lines", "tabbedout", "lines"),
    ])
}

fn exsrch(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    remove(dir, &["q.fa", "q.fq"]);
    generate("exsrch", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let (fmt, o) = (&w[0], &w[1..]);
    let mut vo = Vec::new();
    let mut i = 0;
    while i < o.len() {
        match o[i].as_str() {
            "strand" => vo.extend(strs(&["--strand", "both"])),
            "self" => vo.push("--self".into()),
            n => {
                vo.push(format!("--{n}"));
                if let Some(v) = o.get(i + 1) {
                    vo.push(v.clone());
                }
                i += 1;
            }
        }
        i += 1;
    }
    let mut a = argv(&format!(
        "--search_exact {t}/q.{fmt} --db {t}/db.fa --uc {{O}}/uc --uc_allhits --quiet --notrunclabels --threads 1"
    ));
    a.extend(vo);
    Ok(vec![
        Case::new(format!("fmt={fmt} {}", o.join(" ")), a)
            .errors("abundance")
            .success_only(&["uc"])
            .count("uc", "uc", "uc"),
    ])
}

const HITFMT_FIELDS: &str = "query+target+evalue+id+pctpv+pctgaps+pairs+gaps+qlo+qhi+tlo+thi+pv+ql+tl+qs+ts+\
alnlen+opens+exts+raw+bits+aln+caln+qstrand+tstrand+qrow+trow+qframe+tframe+mism+ids+qcov+tcov+id0+id1+id2+\
id3+id4+qilo+qihi+tilo+tihi+diffs+mid+qseq+tseq+qrowdots+trowdots+qlor+qhir+tlor+thir";

fn hitfmt(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("hitfmt", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let (mask, rowlen, o) = (&w[0], &w[1], &w[2..]);
    let mut a = argv(&format!(
        "--usearch_global {t}/q.fa --db {t}/db.fa --quiet --notrunclabels --threads 1 --qmask {mask} --dbmask {mask}"
    ));
    a.extend(o.iter().cloned());
    a.extend(argv(&format!(
        "--rowlen {rowlen} --uc {{O}}/uc --uc_allhits --alnout {{O}}/aln --blast6out {{O}}/b6 --userout {{O}}/user \
         --userfields {HITFMT_FIELDS} --samout {{O}}/sam --samheader --fastapairs {{O}}/pairs --qsegout {{O}}/qseg \
         --tsegout {{O}}/tseg --fasta_width 0"
    )));
    Ok(vec![
        Case::new(format!("mask={mask} rowlen={rowlen} {}", o.join(" ")), a)
            .success_outputs_only()
            .cover("mask", mask.as_str())
            .count("uc", "uc", "uc"),
    ])
}

fn otutab(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("otutab", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let cmd = if w[0] == "global" {
        "--usearch_global"
    } else {
        "--search_exact"
    };
    let mut a = argv(&format!(
        "{cmd} {t}/q.fa --db {t}/db.fa --quiet --notrunclabels --threads 1 --minseqlength 1"
    ));
    a.extend(w[1..].iter().cloned());
    a.extend(argv(
        "--uc {O}/uc --uc_allhits --otutabout {O}/tsv --mothur_shared_out {O}/mo --biomout {O}/biom",
    ));
    Ok(vec![
        Case::new(format!("{} {}", w[0], w[1..].join(" ")), a)
            .success_outputs_only()
            .cover("mode", w[0].as_str())
            .count("uc", "uc", "uc"),
    ])
}

fn lcatax(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("lcatax", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let mut a = argv(&format!(
        "--usearch_global {t}/q.fa --db {t}/db.fa --quiet --notrunclabels --threads 1 --minseqlength 1"
    ));
    a.extend(w[1..].iter().cloned());
    a.extend(argv(&format!(
        "--lca_cutoff {} --uc {{O}}/uc --uc_allhits --lcaout {{O}}/lca",
        w[0]
    )));
    Ok(vec![
        Case::new(format!("cutoff={} {}", w[0], w[1..].join(" ")), a)
            .success_outputs_only()
            .count("uc", "uc", "uc"),
    ])
}

fn getseq_hits(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("hitfmt", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let (mask, o) = (&w[0], &w[2..]);
    let mut a = argv(&format!(
        "--usearch_global {t}/q.fa --db {t}/db.fa --quiet --notrunclabels --threads 1 --qmask {mask} --dbmask {mask}"
    ));
    a.extend(o.iter().cloned());
    a.extend(argv(
        "--uc {O}/uc --uc_allhits --matched {O}/m --notmatched {O}/n --dbmatched {O}/dm --dbnotmatched {O}/dn \
         --fasta_width 0",
    ));
    Ok(vec![
        Case::new(format!("mask={mask} {}", o.join(" ")), a)
            .success_outputs_only()
            .count("matched", "matched", "records")
            .count("dbmatched", "dbmatched", "records"),
    ])
}

fn allpr(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("allpr", dir, seed, &[])?;
    let t = p(dir);
    let vo = vopts(&params(dir));
    let mut a = argv(&format!(
        "--allpairs_global {t}/in.fa --uc {{O}}/uc --quiet --notrunclabels --qmask none --minseqlength 0 \
         --maxseqlength 100000 --threads 1"
    ));
    a.extend(vo.iter().cloned());
    Ok(vec![
        Case::new(vo.join(" "), a)
            .success_only(&["uc"])
            .count("uc", "uc", "uc"),
    ])
}

fn dbsrch(dir: &Path, seed: u64, env: &Env) -> io::Result<Vec<Case>> {
    generate("dbsrch", dir, seed, &[])?;
    let t = p(dir);
    let mut w = params(dir);
    let mut db = format!("{t}/db.fa");
    if w.first().map(String::as_str) == Some("udb") {
        let udb = format!("{t}/db.udb");
        sh(env, "mkudb", &["-i", &db, "-o", &udb, "--word-len", &w[1]]);
        db = udb;
        w.drain(..2);
    }
    let vo = vopts(&w);
    let mut a = argv(&format!(
        "--usearch_global {t}/q.fa --db {db} --uc {{O}}/uc --uc_allhits --quiet --notrunclabels --qmask soft \
         --dbmask soft --minseqlength 0 --maxseqlength 100000 --threads 1"
    ));
    a.extend(vo.iter().cloned());
    Ok(vec![
        Case::new(vo.join(" "), a)
            .success_only(&["uc"])
            .cover("db", &db[db.len() - 3..])
            .count("uc", "uc", "uc"),
    ])
}

fn cluster_cmd(mode: &str, input: &str) -> String {
    match mode {
        "smallmem" => format!("--cluster_smallmem {input} --usersort"),
        "fast" => format!("--cluster_fast {input}"),
        _ => format!("--cluster_size {input}"),
    }
}

const COMMON: &str = "--quiet --notrunclabels --qmask soft --minseqlength 0 --maxseqlength 100000 --threads 1 \
--fasta_width 0";

/// clust's parameters: the mode, the option words, and whether
/// `nosizein` and `nosizeout` are given.
fn cluster_params(dir: &Path) -> (String, Vec<String>, bool, bool) {
    let w = params(dir);
    let mode = w[0].clone();
    let (mut nosizein, mut nosizeout, mut words) = (false, false, Vec::new());
    for x in &w[1..] {
        match x.as_str() {
            "nosizein" => nosizein = true,
            "nosizeout" => nosizeout = true,
            _ => words.push(x.clone()),
        }
    }
    (mode, words, nosizein, nosizeout)
}

fn sizes(nosizein: bool, nosizeout: bool) -> Vec<String> {
    let mut s = Vec::new();
    if !nosizein {
        s.push("--sizein".to_string());
    }
    if !nosizeout {
        s.push("--sizeout".to_string());
    }
    s
}

fn clust(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("clust", dir, seed, &[])?;
    let (mode, words, nin, nout) = cluster_params(dir);
    let s = sizes(nin, nout);
    let vo = vopts(&words);
    let mut a = argv(&format!(
        "{} --uc {{O}}/uc --centroids {{O}}/c {COMMON}",
        cluster_cmd(&mode, &format!("{}/in.fa", p(dir)))
    ));
    a.extend(s.iter().cloned());
    a.extend(vo.iter().cloned());
    let name = format!("{mode} {}", [vo, s].concat().join(" "));
    Ok(vec![
        Case::new(name, a)
            .success_only(&["uc", "centroids"])
            .cover("mode", mode.as_str())
            .count("uc", "uc", "uc"),
    ])
}

fn unoise(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("unoise", dir, seed, &[])?;
    let w = params(dir);
    let s = sizes(
        w.iter().any(|x| x == "nosizein"),
        w.iter().any(|x| x == "nosizeout"),
    );
    let words: Vec<String> = w
        .into_iter()
        .filter(|x| x != "nosizein" && x != "nosizeout")
        .collect();
    let vo = vopts(&words);
    let mut a = argv(&format!(
        "--cluster_unoise {}/in.fa --uc {{O}}/uc --centroids {{O}}/c --quiet --notrunclabels --qmask soft \
         --minseqlength 1 --maxseqlength 100000 --threads 1 --fasta_width 0",
        p(dir)
    ));
    a.extend(s.iter().cloned());
    a.extend(vo.iter().cloned());
    Ok(vec![
        Case::new([vo, s].concat().join(" "), a)
            .success_only(&["uc", "centroids"])
            .count("uc", "uc", "uc"),
    ])
}

/// The --id, --strand and --wordlength of clust's parameters.
fn id_strand_words(words: &[String]) -> Vec<String> {
    let mut vo = Vec::new();
    let mut i = 0;
    while i < words.len() {
        let w = words[i].as_str();
        if ["id", "strand", "wordlength"].contains(&w) {
            vo.push(format!("--{w}"));
            if let Some(v) = words.get(i + 1) {
                vo.push(v.clone());
            }
            i += 2;
        } else if VFLAGS.contains(&w) {
            i += 1;
        } else {
            i += 2;
        }
    }
    vo
}

fn clspl(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("clust", dir, seed, &[])?;
    let (mode, words, _, _) = cluster_params(dir);
    let vo = id_strand_words(&words);
    let mut a = argv(&format!(
        "{} --uc {{O}}/uc --clusters {{O}}/cluster --sizein {COMMON}",
        cluster_cmd(&mode, &format!("{}/in.fa", p(dir)))
    ));
    a.extend(vo.iter().cloned());
    Ok(vec![
        Case::new(format!("{mode} {}", vo.join(" ")), a)
            .success_only(&["uc", "clusters"])
            .cover("mode", mode.as_str()),
    ])
}

fn clmsa(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("clust", dir, seed, &[])?;
    let (mode, words, nin, nout) = cluster_params(dir);
    let vo = id_strand_words(&words);
    let s = sizes(nin, nout);
    let mut a = argv(&format!(
        "{} --msaout {{O}}/msa --consout {{O}}/cons --profile {{O}}/prof {COMMON}",
        cluster_cmd(&mode, &format!("{}/in.fa", p(dir)))
    ));
    a.extend(s.iter().cloned());
    a.extend(vo.iter().cloned());
    Ok(vec![
        Case::new(format!("{mode} {}", [vo, s].concat().join(" ")), a)
            .success_only(&["msaout", "consout", "profile"])
            .cover("mode", mode.as_str())
            .count("clusters", "consout", "records"),
    ])
}

/// name value pairs as vsearch options.
fn pairs(words: &[String]) -> Vec<String> {
    words
        .iter()
        .enumerate()
        .map(|(k, x)| {
            if k % 2 == 0 {
                format!("--{x}")
            } else {
                x.clone()
            }
        })
        .collect()
}

fn uchime(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("uchime", dir, seed, &[])?;
    let t = p(dir);
    let w = params(dir);
    let mode = w[0].clone();
    let vo = pairs(&w[1..]);
    let db = if mode == "uchime_ref" {
        format!("--db {t}/db.fa --dbmask soft")
    } else {
        String::new()
    };
    let mut a = argv(&format!(
        "--{mode} {t}/in.fa {db} --nonchimeras {{O}}/non --chimeras {{O}}/chim --borderline {{O}}/bord \
         --uchimeout {{O}}/tab --uchimealns {{O}}/aln {COMMON}"
    ));
    a.extend(vo.iter().cloned());
    Ok(vec![
        Case::new(format!("{mode} {}", vo.join(" ")), a)
            .success_only(&[
                "nonchimeras",
                "chimeras",
                "borderline",
                "uchimeout",
                "uchimealns",
            ])
            .cover("mode", mode.as_str())
            .count("Y", "uchimeout", "col 17 Y")
            .count("N", "uchimeout", "col 17 N")
            .count("border", "uchimeout", "col 17 ?"),
    ])
}

fn mchime(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("mchime", dir, seed, &[])?;
    let t = p(dir);
    let vo = pairs(&params(dir));
    let mut a = argv(&format!(
        "--chimeras_denovo {t}/in.fa --nonchimeras {{O}}/non --chimeras {{O}}/chim --tabbedout {{O}}/tab \
         --alnout {{O}}/aln {COMMON}"
    ));
    a.extend(vo.iter().cloned());
    Ok(vec![
        Case::new(vo.join(" "), a)
            .success_only(&["nonchimeras", "chimeras", "tabbedout", "alnout"])
            .count("Y", "tabbedout", "lines"),
    ])
}
