//! Malformed FASTA and FASTQ against the vsearch commands, from the former
//! `input.zsh` and `report-input.zsh`.
//!
//! Each check compares the output, the exit code, the error message
//! (vsearch's "Fatal error:" with rsearch's "Error:") and the warning about
//! stripped characters. Every check runs the translated pipeline. Where
//! that pipeline starts with another tool than the one `input.zsh` fuzzed
//! (conv --fasta, relab, topn, subseq, mkudb, dbornt, allpr), `input.zsh`'s
//! command also runs, as a probe named with "probe".

use std::io;
use std::path::Path;

use super::{Env, Source, argv, generate, p, seeds};
use crate::cases::{Case, Probe};
use crate::translate::Step;

pub fn sources() -> Vec<Source> {
    vec![
        Source {
            name: "fuzz",
            seeds: || seeds(1, 400),
            about: "input checks of every tool on malformed FASTA and FASTQ (input.zsh)",
            make: fuzz,
        },
        Source {
            name: "fuzz-report",
            seeds: || seeds(1, 400),
            about: "the quality reports on malformed FASTQ (report-input.zsh)",
            make: fuzz_report,
        },
    ]
}

fn files(dir: &Path) -> io::Result<()> {
    for (name, text) in [
        (
            "odb.fa",
            ">d1\nACGTACGTNACGTTGCAAGTCCGT\n>d2\nTTGACCAGTAGGCATTAGCAGGTA\n",
        ),
        (
            "oq.fa",
            ">q1\nACGTTGCAAGTC\n>q2\nTACCTGCTAATG\n>q3\nGGGGGGG\n",
        ),
        (
            "sdb.fa",
            ">t1;tax=d:A,g:B\nACGTACGTNACGTTGCAAGTCCGTAACCGGTTACGATCGATCGGCTA\n\
             >t2;tax=d:C\nTTGACCAGTAGGCATTAGCAGGTAACGGATCCATGCAAGGCTTGCA\n",
        ),
        (
            "sq.fa",
            ">q1\nACGTTGCAAGTCCGTAACCGGTTACGATCGATCG\n>q2\nTACCTGCTAATGCCTACTGGTCAAGGCAGGTA\n",
        ),
        ("edb.fa", ">e1\nACGT\n>e2\nacgtn\n>e3\nA\n"),
        ("eq.fa", ">q1\nACGT\n>q2\nA\n>q3\nACGTN\n"),
        (
            "ddb.fa",
            ">d1\nACGTACGTTGCAACGTAGCT\n>d2\nacgtACGTTGCAACGTAGCA\n>d3\nTTGCAACG\n",
        ),
        (
            "dq.fa",
            ">p1\nACGTACGTTGCAACGTAGCT\n>p2\nAGCTACGTTGCAACGTACGT\n>p3\n\n",
        ),
        (
            "udb.fa",
            ">u1;size=9\nACGTACGTTGCAACGTAGCTAGGCTAACGTTAGCAT\n>u2;size=8\nTTGCAACGGTCATTGACCATGGTACCAGTTAACGGA\n",
        ),
    ] {
        std::fs::write(dir.join(name), text)?;
    }
    Ok(())
}

/// `input.zsh`'s rsearch command: steps of (tool, args), each into the next
/// (the first reads its `-i`); the last one's output stands for `out`.
fn chain(steps: Vec<(&'static str, Vec<String>)>, out: &str) -> Probe {
    let out = out.to_string();
    Box::new(move |rdir: &Path| {
        let mut made = Vec::new();
        let mut prev: Option<String> = None;
        for (k, (tool, args)) in steps.iter().enumerate() {
            let dst = p(&rdir.join(format!("p{k}")));
            made.push(Step {
                tool: tool.to_string(),
                args: args.clone(),
                stdin: prev.clone(),
                stdout: Some(dst.clone()),
                system: false,
            });
            prev = Some(dst);
        }
        (made, vec![(out.clone(), prev.unwrap())])
    })
}

fn s(words: &[&str]) -> Vec<String> {
    words.iter().map(|w| w.to_string()).collect()
}

fn fuzz(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("fuzz", dir, seed, &[])?;
    files(dir)?;
    let t = p(dir);
    let q = "--quiet --notrunclabels";
    let mut out: Vec<Case> = Vec::new();
    // a check: the translated pipeline, and input.zsh's command as a probe
    let chk = |out: &mut Vec<Case>,
               name: &str,
               v: &str,
               probe: Option<Probe>,
               ok_on_success: bool,
               errors: &str| {
        let a = argv(&v.replace("{T}", &t));
        let outputs_on_error = !ok_on_success;
        let mut c = Case::new(name, a.clone())
            .errors(errors)
            .warn(&["stripped"]);
        // where the pipeline streams into a tool other than the one fuzzed,
        // what it writes before an error is the probe's to check
        if !(outputs_on_error && probe.is_none()) {
            c = c.success_outputs_only();
        }
        out.push(c);
        if let Some(probe) = probe {
            let mut c = Case::new(format!("{name} probe"), a)
                .errors(errors)
                .warn(&["stripped"])
                .probe(probe);
            if !outputs_on_error {
                c = c.success_outputs_only();
            }
            out.push(c);
        }
    };
    for f in ["fq", "fa"] {
        let file = format!("{t}/f.{f}");
        let (o, oo) = if f == "fq" {
            ("--fastqout", "fastqout")
        } else {
            ("--fasta_width 0 --fastaout", "fastaout")
        };
        let cmd = if f == "fq" {
            "--fastq_convert"
        } else {
            "--fasta2fastq"
        };
        let mut relab = vec![("relab", s(&["-i", &file]))];
        if f == "fa" {
            relab.push(("conv", vec![]));
        }
        // conv's translation, with relab's input checks as the probe
        chk(
            &mut out,
            &format!("conv-{f}"),
            &format!("{cmd} {file} --quiet --fastqout {{O}}/out"),
            Some(chain(relab, "fastqout")),
            false,
            "exact",
        );
        out.last_mut().unwrap().name = format!("relab-{f} probe");
        chk(
            &mut out,
            &format!("conv-fasta-{f}"),
            &format!(
                "--fastx_subsample {file} --quiet --sample_pct 100 --minseqlength 0 --notrunclabels --fastaout {{O}}/out"
            ),
            Some(chain(
                vec![("conv", s(&["-i", &file, "--fasta"]))],
                "fastaout",
            )),
            true,
            "exact",
        );
        chk(
            &mut out,
            &format!("topn-{f}"),
            &format!(
                "--sortbysize {file} {q} --minseqlength 0 --fasta_width 0 --topn 1000000 --output {{O}}/out"
            ),
            None,
            false,
            "exact",
        );
        out.push(
            Case::new(
                format!("topn-{f} probe"),
                argv(&format!(
                    "--fastx_subsample {file} {q} --sample_pct 100 --minseqlength 0 {o} {{O}}/out"
                )),
            )
            .errors("exact")
            .warn(&["stripped"])
            .success_outputs_only()
            .probe(chain(vec![("topn", s(&["-i", &file, "1000000"]))], oo)),
        );
        chk(
            &mut out,
            &format!("derep-{f}"),
            &format!(
                "--fastx_uniques {file} --quiet --minseqlength 0 --notrunclabels --sizein --sizeout {o} {{O}}/out"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("sqfilt-{f}"),
            &format!("--fastx_filter {file} {q} {o} {{O}}/out"),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("revcomp-{f}"),
            &format!("--fastx_revcomp {file} --quiet {o} {{O}}/out"),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("syncpa-fwd-{f}"),
            &format!("--fastx_syncpairs {file} --reverse /dev/null --quiet {o}_orphans {{O}}/out"),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("syncpa-rev-{f}"),
            &format!(
                "--fastx_syncpairs /dev/null --reverse {file} --quiet {o}_orphans_rev {{O}}/out"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("lcmask-{f}"),
            &format!("--fastx_mask {file} {q} --minseqlength 0 {o} {{O}}/out"),
            None,
            true,
            "exact",
        );
        chk(
            &mut out,
            &format!("getseq-{f}"),
            &format!("--fastx_getseq {file} {q} --label_substr_match --label a {o} {{O}}/out"),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("subseq-{f}"),
            &format!(
                "--fastx_getsubseq {file} {q} --label_substr_match --label '' --subseq_start 2 --subseq_end 9 {o} {{O}}/out"
            ),
            Some(chain(
                vec![
                    ("subseq", s(&["--start", "2", "--end", "9", "-i", &file])),
                    ("getseq", s(&["--substr", "--label", ""])),
                ],
                oo,
            )),
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("shufl-{f}"),
            &format!(
                "--shuffle {file} {q} --minseqlength 0 --fasta_width 0 --randseed 7 --output {{O}}/out"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("subsmp-{f}"),
            &format!(
                "--fastx_subsample {file} {q} --minseqlength 0 --randseed 3 --sample_pct 60 {o} {{O}}/out"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("scrmbl-{f}"),
            &format!("--scramble {file} {q} --randseed 5 --scramble_kmer 3 {o} {{O}}/out"),
            None,
            false,
            "exact",
        );
        // mkudb reads everything before writing. As vsearch upper-cases,
        // the probe compares the UDB file of the input through lcmask
        // --clear
        let mk_file = file.clone();
        let mk: Probe = Box::new(move |rdir: &Path| {
            let (p0, p1, p2) = (
                p(&rdir.join("p0")),
                p(&rdir.join("p1")),
                p(&rdir.join("p2")),
            );
            let st = |tool: &str, args: Vec<String>, stdin: Option<String>, stdout: String| Step {
                tool: tool.into(),
                args,
                stdin,
                stdout: Some(stdout),
                system: false,
            };
            (
                vec![
                    st("mkudb", s(&["-i", &mk_file, "--word-len", "3"]), None, p0),
                    st("lcmask", s(&["--clear", "-i", &mk_file]), None, p1.clone()),
                    st("mkudb", s(&["--word-len", "3"]), Some(p1), p2.clone()),
                ],
                vec![("output".to_string(), p2)],
            )
        });
        chk(
            &mut out,
            &format!("mkudb-{f}"),
            &format!(
                "--makeudb_usearch {file} {q} --minseqlength 0 --dbmask none --wordlength 3 --output {{O}}/out"
            ),
            Some(mk),
            true,
            "mkudb-exact",
        );
        chk(
            &mut out,
            &format!("dbornt-{f}"),
            &format!("--orient {file} --db {{T}}/odb.fa --wordlength 3 {q} {o} {{O}}/out"),
            Some(chain(
                vec![(
                    "dbornt",
                    s(&[
                        "-i",
                        &file,
                        "--db",
                        &format!("{t}/odb.fa"),
                        "--word-len",
                        "3",
                    ]),
                )],
                oo,
            )),
            false,
            "exact",
        );
        let (oq, db_file) = (format!("{t}/oq.fa"), file.clone());
        let db_probe: Probe = Box::new(move |rdir: &Path| {
            let tab = p(&rdir.join("tab"));
            (
                vec![Step {
                    tool: "dbornt".into(),
                    args: s(&[
                        "-i",
                        &oq,
                        "--db",
                        &db_file,
                        "--word-len",
                        "3",
                        "--tabbedout",
                        &tab,
                    ]),
                    stdin: None,
                    stdout: Some(p(&rdir.join("p0"))),
                    system: false,
                }],
                vec![("tabbedout".to_string(), tab)],
            )
        });
        chk(
            &mut out,
            &format!("dbornt-db-{f}"),
            &format!(
                "--orient {{T}}/oq.fa --db {file} --dbmask soft --wordlength 3 {q} --minseqlength 1 --tabbedout {{O}}/out"
            ),
            Some(db_probe),
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("sintax-{f}"),
            &format!(
                "--sintax {file} --db {{T}}/sdb.fa --wordlength 3 --randseed 3 --threads 1 {q} --minseqlength 0 \
                 --tabbedout {{O}}/out"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("sintax-db-{f}"),
            &format!(
                "--sintax {{T}}/sq.fa --db {file} --dbmask soft --wordlength 3 --randseed 3 --threads 1 {q} \
                 --minseqlength 0 --sintax_random --tabbedout {{O}}/out"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("exsrch-{f}"),
            &format!(
                "--search_exact {file} --db {{T}}/edb.fa --uc {{O}}/out --uc_allhits --strand both {q} --threads 1"
            ),
            None,
            true,
            "exact",
        );
        chk(
            &mut out,
            &format!("exsrch-db-{f}"),
            &format!(
                "--search_exact {{T}}/eq.fa --db {file} --uc {{O}}/out --uc_allhits {q} --threads 1 --minseqlength 0"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("allpr-{f}"),
            &format!(
                "--allpairs_global {file} --acceptall {q} --qmask none --minseqlength 0 --threads 1 --uc {{O}}/out"
            ),
            Some(chain(
                vec![("allpr", s(&["-i", &file, "--accept-all"]))],
                "uc",
            )),
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("dbsrch-{f}"),
            &format!(
                "--usearch_global {file} --db {{T}}/ddb.fa --id 0.5 --wordlength 3 --maxaccepts 0 --maxrejects 0 \
                 --strand both --uc {{O}}/out --uc_allhits --qmask soft --dbmask soft {q} --minseqlength 0 --threads 1"
            ),
            None,
            true,
            "exact",
        );
        chk(
            &mut out,
            &format!("dbsrch-db-{f}"),
            &format!(
                "--usearch_global {{T}}/dq.fa --db {file} --id 0.5 --wordlength 3 --maxaccepts 4 --uc {{O}}/out \
                 --uc_allhits --qmask soft --dbmask soft {q} --minseqlength 0 --maxseqlength 1000000000 --threads 1"
            ),
            None,
            false,
            "exact",
        );
        let big = "--minseqlength 0 --maxseqlength 1000000000 --threads 1";
        chk(
            &mut out,
            &format!("clust-{f}"),
            &format!(
                "--cluster_smallmem {file} --usersort --id 0.8 --wordlength 3 --strand both --uc {{O}}/out --sizein \
                 {q} --qmask soft {big}"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("unoise-{f}"),
            &format!(
                "--cluster_unoise {file} --minsize 1 --wordlength 3 --uc {{O}}/out --sizein {q} --qmask soft {big}"
            ),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("uchime-{f}"),
            &format!("--uchime_denovo {file} --uchimeout {{O}}/out {q} --qmask soft {big}"),
            None,
            false,
            "exact",
        );
        chk(
            &mut out,
            &format!("mchime-{f}"),
            &format!("--chimeras_denovo {file} --tabbedout {{O}}/out {q} --qmask soft {big}"),
            None,
            false,
            "exact",
        );
    }
    let f = format!("{t}/f");
    chk(
        &mut out,
        "rerep-fa",
        &format!("--rereplicate {f}.fa {q} --sizeout --fasta_width 0 --output {{O}}/out"),
        None,
        false,
        "exact",
    );
    chk(
        &mut out,
        "joinpa-fq",
        &format!("--fastq_join {f}.fq --reverse {f}.fq --quiet --fastqout {{O}}/out"),
        None,
        false,
        "exact",
    );
    chk(
        &mut out,
        "mergpa-fq",
        &format!("--fastq_mergepairs {f}.fq --reverse {f}.fq --quiet --fastqout {{O}}/out"),
        None,
        true,
        "exact",
    );
    chk(
        &mut out,
        "redig-fa",
        &format!(
            "--cut {f}.fa --cut_pattern G^AATT_C --quiet --fasta_width 0 --fastaout {{O}}/out"
        ),
        None,
        false,
        "exact",
    );
    chk(
        &mut out,
        "sortby-fq",
        &format!("--sortbylength {f}.fq {q} --minseqlength 0 --fasta_width 0 --output {{O}}/out"),
        None,
        false,
        "exact",
    );
    chk(
        &mut out,
        "sortby-fa",
        &format!("--sortbysize {f}.fa {q} --minseqlength 0 --fasta_width 0 --output {{O}}/out"),
        None,
        false,
        "exact",
    );
    chk(
        &mut out,
        "uchime-ref-fa",
        &format!(
            "--uchime_ref {f}.fa --db {{T}}/udb.fa --uchimeout {{O}}/out {q} --qmask soft --dbmask soft \
             --minseqlength 0 --threads 1"
        ),
        None,
        true,
        "exact",
    );
    Ok(out)
}

fn fuzz_report(dir: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    generate("fuzz", dir, seed, &[])?;
    let t = p(dir);
    Ok([
        ("fastq_chars", "--log"),
        ("fastq_stats", "--log"),
        ("fastq_eestats", "--output"),
        ("fastq_eestats2", "--output"),
    ]
    .iter()
    .map(|(cmd, o)| {
        Case::new(
            *cmd,
            argv(&format!("--{cmd} {t}/f.fq --quiet {o} {{O}}/out")),
        )
        .errors("exact")
        .success_only(&[o.trim_start_matches('-')])
    })
    .collect())
}
