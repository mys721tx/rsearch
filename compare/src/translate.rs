//! Translates a vsearch command line into the rsearch pipeline that gives
//! the same results.
//!
//! The rules come from the pipelines of `README.md`, the tool pages in
//! `docs/` and the option map in `docs/vsearch.md`, each written once here.
//! [`translate`] returns a [`Pipeline`]. It holds the rsearch commands in
//! the order they run ([`Step`]s, each a tool, its arguments and the files
//! on its standard input and output). For each output of the vsearch
//! command, it also holds the rsearch file that stands for it and the
//! comparator that compares them ([`Output`]s). Temporary files and
//! rsearch's outputs go into the working directory. Steps run one after the
//! other, each reading the files the earlier ones wrote, which is a shell
//! pipeline without the concurrency.
//!
//! An option it cannot translate is a [`TranslationError`] naming it. A
//! behaviour that is vsearch's default rather than an option, and that the
//! pipeline does not reproduce (a length filter sqfilt cannot express, labels
//! cut at spaces without `--notrunclabels`, and the order threaded vsearch
//! writes queries in), is listed in `caveats` instead.

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::process::{Command, Stdio};

pub use crate::options::TranslationError;
use crate::options::{Args, err};

/// One command of a pipeline.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Step {
    pub tool: String,
    pub args: Vec<String>,
    pub stdin: Option<String>,
    pub stdout: Option<String>,
    /// A program on `PATH` (gzip, bzip2), not one of rsearch's tools.
    pub system: bool,
}

impl Step {
    /// The program and its arguments, rsearch's tools taken from `bindir`.
    pub fn argv(&self, bindir: &Path) -> Vec<String> {
        let program = if self.system {
            self.tool.clone()
        } else {
            bindir.join(&self.tool).to_string_lossy().into_owned()
        };
        std::iter::once(program)
            .chain(self.args.iter().cloned())
            .collect()
    }

    /// The files this step reads: its standard input, `-i` and `--reverse`.
    pub fn inputs(&self) -> Vec<String> {
        let mut out: Vec<String> = self.stdin.iter().cloned().collect();
        for k in 0..self.args.len().saturating_sub(1) {
            if self.args[k] == "-i" || self.args[k] == "--reverse" {
                out.push(self.args[k + 1].clone());
            }
        }
        out
    }
}

/// A vsearch output, and the rsearch file that stands for it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Output {
    /// The vsearch option, without dashes (`fastaout`).
    pub option: String,
    /// Its file (or prefix, for `--clusters`) in vsearch's argv.
    pub vsearch: String,
    /// The rsearch file (or prefix) that stands for it.
    pub rsearch: String,
    /// The comparator that compares them (see `comparators`).
    pub comparator: String,
}

/// The translation of a vsearch command line.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Pipeline {
    pub command: String,
    pub argv: Vec<String>,
    pub steps: Vec<Step>,
    pub outputs: Vec<Output>,
    pub caveats: Vec<String>,
}

impl Pipeline {
    /// runs the steps of rsearch's threaded tools on `threads` threads.
    pub fn set_threads(&mut self, threads: usize) {
        for step in &mut self.steps {
            if step.system || !THREADED_TOOLS.contains(&step.tool.as_str()) {
                continue;
            }
            let mut args = Vec::with_capacity(step.args.len() + 2);
            let mut given = std::mem::take(&mut step.args).into_iter();
            while let Some(a) = given.next() {
                if a == "-t" || a == "--threads" {
                    given.next();
                } else {
                    args.push(a);
                }
            }
            args.extend(["-t".to_string(), threads.to_string()]);
            step.args = args;
        }
    }

    /// The steps as shell pipelines, a line each.
    pub fn render(&self) -> String {
        let mut uses: HashMap<&str, usize> = HashMap::new();
        let inputs: Vec<Vec<String>> = self.steps.iter().map(Step::inputs).collect();
        for list in &inputs {
            for f in list {
                *uses.entry(f.as_str()).or_default() += 1;
            }
        }
        let mut lines = Vec::new();
        let mut pending: Option<(String, Option<String>)> = None;
        for s in &self.steps {
            let words: Vec<String> = std::iter::once(s.tool.clone())
                .chain(s.args.iter().map(|a| shell_quote(a)))
                .collect();
            let joins = match (&pending, &s.stdin) {
                (Some((_, out)), Some(stdin)) => {
                    out.as_deref() == Some(stdin.as_str())
                        && uses.get(stdin.as_str()) == Some(&1)
                        && is_tmp(stdin)
                }
                _ => false,
            };
            let text = if joins {
                format!("{} | {}", pending.take().unwrap().0, words.join(" "))
            } else {
                if let Some(p) = pending.take() {
                    lines.push(close(p));
                }
                let mut text = words.join(" ");
                if let Some(stdin) = &s.stdin {
                    text.push_str(" < ");
                    text.push_str(stdin);
                }
                text
            };
            pending = Some((text, s.stdout.clone()));
        }
        if let Some(p) = pending {
            lines.push(close(p));
        }
        lines.join("\n")
    }

    /// The pipeline as JSON, the shape `rscompare translate --json` prints.
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "command": self.command,
            "steps": self.steps.iter().map(|s| serde_json::json!({
                "tool": s.tool, "args": s.args, "stdin": s.stdin, "stdout": s.stdout,
            })).collect::<Vec<_>>(),
            "outputs": self.outputs.iter().map(|o| serde_json::json!({
                "option": o.option, "vsearch": o.vsearch, "rsearch": o.rsearch,
                "comparator": o.comparator,
            })).collect::<Vec<_>>(),
            "caveats": self.caveats,
        })
    }
}

fn close((text, out): (String, Option<String>)) -> String {
    match out {
        Some(out) => format!("{text} > {out}"),
        None => text,
    }
}

/// A word as the shell reads it.
pub fn shell_quote(a: &str) -> String {
    if a.is_empty() || a.chars().any(|c| " \t'\"$*;|&<>()^\\".contains(c)) {
        format!("'{}'", a.replace('\'', "'\\''"))
    } else {
        a.to_string()
    }
}

/// Whether a path is one of the translator's temporary files.
pub fn is_tmp(path: &str) -> bool {
    Path::new(path)
        .file_name()
        .is_some_and(|n| n.to_string_lossy().starts_with("tmp"))
}

/// The compression of a file by its magic number, as vsearch detects it:
/// `gzip`, `bzip2`, or none.
pub fn compression(path: &str) -> Option<&'static str> {
    if path.is_empty() || path == "-" {
        return None;
    }
    let mut head = [0u8; 3];
    let n = read_head(path, &mut head)?;
    if n >= 2 && head[..2] == [0x1f, 0x8b] {
        Some("gzip")
    } else if n == 3 && &head == b"BZh" {
        Some("bzip2")
    } else {
        None
    }
}

fn read_head(path: &str, buf: &mut [u8]) -> Option<usize> {
    let mut f = std::fs::File::open(path).ok()?;
    let mut n = 0;
    while n < buf.len() {
        match f.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(_) => return None,
        }
    }
    Some(n)
}

/// The format of a file by its first bytes, inside its compression:
/// `fasta`, `fastq`, `udb`, `sff`, or none.
pub fn sniff(path: &str) -> Option<&'static str> {
    if path.is_empty() || path == "-" {
        return None;
    }
    let mut head = [0u8; 4];
    let n = match compression(path) {
        None => read_head(path, &mut head)?,
        Some(kind) => {
            let mut child = Command::new(kind)
                .args(["-dc", path])
                .stdout(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()
                .ok()?;
            let mut out = child.stdout.take()?;
            let mut n = 0;
            while n < head.len() {
                match out.read(&mut head[n..]) {
                    Ok(0) | Err(_) => break,
                    Ok(k) => n += k,
                }
            }
            drop(out);
            let _ = child.kill();
            let _ = child.wait();
            n
        }
    };
    match &head[..n] {
        [b'>', ..] => Some("fasta"),
        [b'@', ..] => Some("fastq"),
        b"FBDU" => Some("udb"),
        b".sff" => Some("sff"),
        _ => None,
    }
}

/// Acceptance criteria and scores: vsearch name and rsearch option.
const CRITERIA: &[(&str, &str)] = &[
    ("id", "--id"),
    ("iddef", "--iddef"),
    ("weak_id", "--weak-id"),
    ("maxdiffs", "--max-diffs"),
    ("maxsubs", "--max-subs"),
    ("maxgaps", "--max-gaps"),
    ("mincols", "--min-cols"),
    ("query_cov", "--query-cov"),
    ("target_cov", "--target-cov"),
    ("maxid", "--max-id"),
    ("mid", "--mid"),
    ("idprefix", "--id-prefix"),
    ("idsuffix", "--id-suffix"),
    ("maxqt", "--max-qt"),
    ("minqt", "--min-qt"),
    ("minsl", "--min-sl"),
    ("maxsl", "--max-sl"),
    ("maxqsize", "--max-qsize"),
    ("mintsize", "--min-tsize"),
    ("minsizeratio", "--min-size-ratio"),
    ("maxsizeratio", "--max-size-ratio"),
    ("match", "--match"),
    ("mismatch", "--mismatch"),
    ("gapopen", "--gapopen"),
    ("gapext", "--gapext"),
];
const CRITERIA_FLAGS: &[(&str, &str)] = &[
    ("leftjust", "--left-just"),
    ("rightjust", "--right-just"),
    ("self", "--self"),
    ("selfid", "--self-id"),
    ("n_mismatch", "--n-mismatch"),
];
const WORDS: &[(&str, &str)] = &[
    ("strand", "--strand"),
    ("maxaccepts", "--maxaccepts"),
    ("maxrejects", "--maxrejects"),
    ("wordlength", "--word-len"),
    ("minwordmatches", "--min-word-matches"),
];
/// Options vsearch accepts and ignores.
const IGNORED: &[&str] = &[
    "band", "hspw", "minhsp", "xdrop_nw", "fulldp", "slots", "pattern",
];
/// hitfmt takes the options the search scored with. The UC hits do not hold
/// them, and it recomputes the scores and counts from the CIGAR with them.
const HITFMT_SCORING: &[(&str, &str)] = &[
    ("iddef", "--iddef"),
    ("match", "--match"),
    ("mismatch", "--mismatch"),
    ("gapopen", "--gapopen"),
    ("gapext", "--gapext"),
];
const LABELS: &[(&str, &str)] = &[
    ("relabel", "--relabel"),
    ("sample", "--sample"),
    ("label_suffix", "--label-suffix"),
];
const LABEL_FLAGS: &[(&str, &str)] = &[
    ("relabel_keep", "--relabel-keep"),
    ("relabel_md5", "--relabel-md5"),
    ("relabel_self", "--relabel-self"),
    ("relabel_sha1", "--relabel-sha1"),
    ("sizeout", "--sizeout"),
    ("xsize", "--xsize"),
    ("lengthout", "--lengthout"),
    ("xlength", "--xlength"),
    ("xee", "--xee"),
    ("eeout", "--ee"),
    ("fastq_eeout", "--ee"),
];
/// Commands that write their per-query outputs in the order their threads
/// finish when they run on several.
const UNORDERED: &[&str] = &[
    "usearch_global",
    "search_exact",
    "allpairs_global",
    "uchime_ref",
    "sintax",
];
/// rsearch's tools that take `-t/--threads`.
pub const THREADED_TOOLS: &[&str] = &[
    "allpr", "dbornt", "dbsrch", "lcmask", "mergpa", "sintax", "uchime",
];
/// Commands that keep whole headers without `--notrunclabels`, or read none.
pub const WHOLE_HEADERS: &[&str] = &[
    "fastq_convert",
    "fasta2fastq",
    "fastx_filter",
    "fastq_filter",
    "fastx_revcomp",
    "sintax",
    "fastq_chars",
    "fastq_stats",
    "fastq_eestats",
    "fastq_eestats2",
    "udbinfo",
    "udbstats",
    "udb2fasta",
    "sff_convert",
    "fastq_join",
    "fastq_mergepairs",
    "fastx_syncpairs",
    "cut",
    "scramble",
];

type R<T> = Result<T, TranslationError>;
type Cmds = Vec<(&'static str, Vec<String>)>;

fn s(words: &[&str]) -> Vec<String> {
    words.iter().map(|w| w.to_string()).collect()
}

/// A number as Python's `str` of what `int` or `float` makes of it, which
/// the translator's messages and arguments show.
fn py_num(text: &str) -> String {
    if text.contains(['.', 'e', 'E']) {
        match text.trim().parse::<f64>() {
            Ok(v) => py_float(v),
            Err(_) => text.to_string(),
        }
    } else {
        match text.trim().parse::<i64>() {
            Ok(v) => v.to_string(),
            Err(_) => text.to_string(),
        }
    }
}

/// Python's `repr` of a float.
pub fn py_float(v: f64) -> String {
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let e = format!("{v:e}");
    let (mantissa, exp) = e.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    if (-4..16).contains(&exp) {
        let d = format!("{v}");
        if d.contains('.') { d } else { format!("{d}.0") }
    } else {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exp.abs())
    }
}

/// A number given (or defaulted), with the text Python would show.
#[derive(Clone)]
struct Num {
    v: f64,
    text: String,
}

struct Translator {
    a: Args,
    cmd: String,
    inp: String,
    work: String,
    p: Pipeline,
    n: usize,
    width: Num,
    formats: HashMap<String, Option<&'static str>>,
}

impl Translator {
    fn fmt(&self, path: &str) -> Option<&'static str> {
        match self.formats.get(path) {
            Some(f) => *f,
            None => sniff(path),
        }
    }

    /// A compressed input through gzip -dc or bzip2 -dc. vsearch knows a
    /// compressed file by its magic number; --gzip_decompress and
    /// --bzip2_decompress matter only for a pipe, and an uncompressed file
    /// is read as it is with them.
    fn decompressed(&mut self, path: String) -> String {
        if path.is_empty() || path == "-" {
            return path;
        }
        let Some(kind) = compression(&path) else {
            return path;
        };
        let fmt = sniff(&path);
        let out = self.tmp(fmt.unwrap_or(kind));
        self.p.steps.push(Step {
            tool: kind.to_string(),
            args: vec!["-dc".into(), path],
            stdin: None,
            stdout: Some(out.clone()),
            system: true,
        });
        self.formats.insert(out.clone(), fmt);
        out
    }

    // building blocks

    fn caveat(&mut self, text: String) {
        if !self.p.caveats.contains(&text) {
            self.p.caveats.push(text);
        }
    }

    fn tmp(&mut self, name: &str) -> String {
        self.n += 1;
        self.path(&format!("tmp{}.{name}", self.n))
    }

    fn path(&self, name: &str) -> String {
        Path::new(&self.work)
            .join(name)
            .to_string_lossy()
            .into_owned()
    }

    /// Adds a step; returns the file its standard output goes to. A
    /// temporary file given first as -i is read on standard input instead,
    /// so that it streams as in a shell pipeline.
    fn run(
        &mut self,
        tool: &str,
        mut args: Vec<String>,
        mut stdin: Option<String>,
        stdout: Option<String>,
        name: &str,
    ) -> String {
        if stdin.is_none() && args.len() >= 2 && args[0] == "-i" && is_tmp(&args[1]) {
            stdin = Some(args[1].clone());
            args.drain(..2);
        }
        let stdout = stdout.unwrap_or_else(|| self.tmp(name));
        self.p.steps.push(Step {
            tool: tool.to_string(),
            args,
            stdin,
            stdout: Some(stdout.clone()),
            system: false,
        });
        stdout
    }

    /// Runs the commands one into the next, from src.
    fn chain(&mut self, mut src: String, cmds: Cmds) -> String {
        for (tool, args) in cmds {
            src = self.step(tool, &src, args, None, tool);
        }
        src
    }

    /// A step reading src: on its standard input when src is a temporary
    /// file, else as -i.
    fn step(
        &mut self,
        tool: &str,
        src: &str,
        args: Vec<String>,
        stdout: Option<String>,
        name: &str,
    ) -> String {
        if is_tmp(src) {
            self.run(tool, args, Some(src.to_string()), stdout, name)
        } else {
            let mut all = vec!["-i".to_string(), src.to_string()];
            all.extend(args);
            self.run(tool, all, None, stdout, name)
        }
    }

    fn output(&mut self, option: &str, rpath: &str, comparator: &str) {
        let vsearch = self.a.take(option).unwrap_or_default();
        self.p.outputs.push(Output {
            option: option.to_string(),
            vsearch,
            rsearch: rpath.to_string(),
            comparator: comparator.to_string(),
        });
    }

    // outputs of records

    /// vsearch's FASTA output: conv --fasta at --fasta_width, unless the
    /// stream is FASTA already and unwrapped.
    fn fasta(&mut self, option: &str, src: &str, fmt: Option<&str>) {
        if !self.a.has(option) {
            return;
        }
        if fmt == Some("fasta") && self.width.v == 0.0 {
            self.output(option, src, "bytes");
            return;
        }
        let mut args = s(&["--fasta"]);
        if self.width.v != 80.0 {
            args.extend(["--fasta-width".to_string(), self.width.text.clone()]);
        }
        let out = if is_tmp(src) {
            self.run("conv", args, Some(src.to_string()), None, option)
        } else {
            let mut all = vec!["-i".to_string(), src.to_string()];
            all.extend(args);
            self.run("conv", all, None, None, option)
        };
        self.output(option, &out, "bytes");
    }

    fn fastq(&mut self, option: &str, src: &str) {
        if self.a.has(option) {
            self.output(option, src, "bytes");
        }
    }

    fn records(&mut self, src: &str, fmt: Option<&str>, fa: Option<&str>, fq: Option<&str>) {
        if let Some(fa) = fa {
            self.fasta(fa, src, fmt);
        }
        if let Some(fq) = fq {
            self.fastq(fq, src);
        }
    }

    /// Outputs rsearch writes unwrapped, which conv cannot wrap.
    fn unwrapped_only(&self, option: &str) -> R<()> {
        if self.a.has(option) && self.width.v != 0.0 {
            return err(format!(
                "--{option} with --fasta_width {}: rsearch writes it unwrapped (give \
                 --fasta_width 0)",
                self.width.text
            ));
        }
        Ok(())
    }

    // options shared by many commands

    /// The value of an option, or a default, as Python's `take(v, default)`.
    fn take_or(&mut self, name: &str, default: &str) -> String {
        self.a.take(name).unwrap_or_else(|| default.to_string())
    }

    fn take_num(&mut self, name: &str, default: Option<f64>) -> R<Option<Num>> {
        match self.a.take(name) {
            None => Ok(default.map(|v| Num {
                v,
                text: py_num(&format!("{v}")),
            })),
            Some(t) => match t.trim().parse::<f64>() {
                Ok(v) => Ok(Some(Num {
                    v,
                    text: py_num(&t),
                })),
                Err(_) => err(format!("--{name} {t}: not a number")),
            },
        }
    }

    /// The options mapped one to one, taken in table order.
    fn mapped(&mut self, table: &[(&str, &str)]) -> Vec<String> {
        let mut out = Vec::new();
        for (v, r) in table {
            if let Some(t) = self.a.take(v) {
                out.extend([r.to_string(), t]);
            }
        }
        out
    }

    fn qual_in(&mut self) -> Vec<String> {
        self.mapped(&[
            ("fastq_ascii", "--ascii"),
            ("fastq_qmin", "--qmin"),
            ("fastq_qmax", "--qmax"),
        ])
    }

    fn qual_out(&mut self) -> Vec<String> {
        self.mapped(&[
            ("fastq_asciiout", "--ascii-out"),
            ("fastq_qminout", "--qmin-out"),
            ("fastq_qmaxout", "--qmax-out"),
        ])
    }

    fn opt_or(&self, name: &str, default: &str) -> String {
        self.a.peek(name).unwrap_or(default).to_string()
    }

    fn label_args(&mut self, offset: Option<&str>, skip: &[&str]) -> Vec<String> {
        let mut out = Vec::new();
        for (v, r) in LABELS {
            if !skip.contains(v) && self.a.has(v) {
                let t = self.a.take(v).unwrap();
                out.extend([r.to_string(), t]);
            }
        }
        for (v, r) in LABEL_FLAGS {
            if !skip.contains(v) && self.a.flag(v) && !out.iter().any(|o| o == r) {
                out.push(r.to_string());
            }
        }
        if out.iter().any(|o| o == "--ee")
            && let Some(offset) = offset
            && offset != "33"
        {
            out.extend(["--ascii".to_string(), offset.to_string()]);
        }
        out
    }

    /// relab for the label options, `skip` left to an earlier step. relab
    /// drops `;size=` and writes it again only with `--sizeout`. A
    /// `--sizeout` left to the earlier step therefore goes to relab too.
    fn labels(&mut self, src: &str, offset: Option<&str>, skip: &[&str]) -> String {
        let mut args = self.label_args(offset, skip);
        if !args.is_empty()
            && skip.contains(&"sizeout")
            && self.a.flag("sizeout")
            && !args.iter().any(|a| a == "--sizeout")
        {
            args.push("--sizeout".into());
        }
        if args.is_empty() {
            src.to_string()
        } else {
            self.run("relab", args, Some(src.to_string()), None, "relab")
        }
    }

    fn no_labels(&self, why: &str, skip: &[&str]) -> R<()> {
        for (v, _) in LABELS.iter().chain(LABEL_FLAGS) {
            if !skip.contains(v) && self.a.has(v) && !self.a.used.contains(*v) {
                return err(format!("--{v} with --{}: {why}", self.cmd));
            }
        }
        Ok(())
    }

    fn lengths(&mut self, min_default: f64) -> R<(Num, Num)> {
        let mn = self.take_num("minseqlength", Some(min_default))?.unwrap();
        let mx = self.take_num("maxseqlength", Some(50000.0))?.unwrap();
        Ok((mn, mx))
    }

    fn no_max_caveat(&mut self, mn: &Num, mx: &Num, what: &str) {
        self.caveat(format!(
            "--maxseqlength {} not applied to the {what}: with --minseqlength {} the empty \
             sequences are kept, which sqfilt cannot do",
            mx.text, mn.text
        ));
    }

    /// sqfilt's arguments for vsearch's length filters (and any other
    /// filter of the same step), or None for no sqfilt.
    fn length_filter(
        &mut self,
        min_default: f64,
        extra: Vec<String>,
        what: &str,
    ) -> R<Option<Vec<String>>> {
        let (mn, mx) = self.lengths(min_default)?;
        if mn.v <= 0.0 {
            if !extra.is_empty() {
                return err(format!(
                    "--minseqlength {} with {}: sqfilt cannot keep empty sequences",
                    mn.text,
                    extra.join(" ")
                ));
            }
            self.no_max_caveat(&mn, &mx, what);
            return Ok(None);
        }
        let mut args = if mn.v == 1.0 {
            vec![]
        } else {
            vec!["--min-len".into(), mn.text]
        };
        args.extend(["--max-len".to_string(), mx.text]);
        args.extend(extra);
        Ok(Some(args))
    }

    fn filtered(&mut self, src: &str, min_default: f64, extra: Vec<String>) -> R<String> {
        match self.length_filter(min_default, extra, "input")? {
            None => Ok(src.to_string()),
            Some(args) => {
                let mut all = vec!["-i".to_string(), src.to_string()];
                all.extend(args);
                Ok(self.run("sqfilt", all, None, None, "sqfilt"))
            }
        }
    }

    /// lcmask runs for --qmask or --dbmask (and --hardmask). The analysis
    /// tools take no masking options and treat lower case as masked.
    fn masking(&mut self, option: &str, style: &str) -> R<Cmds> {
        let mode = self.take_or(option, "dust");
        let hard = self.a.flag("hardmask");
        if !["dust", "soft", "none"].contains(&mode.as_str()) {
            return err(format!("--{option} {mode}: not dust, soft or none"));
        }
        if style == "udb" {
            // --makeudb_usearch upper-cases first
            return Ok(if mode == "dust" {
                if hard {
                    vec![("lcmask", s(&["--clear"])), ("lcmask", s(&["--hard"]))]
                } else {
                    vec![("lcmask", vec![])]
                }
            } else {
                vec![("lcmask", s(&["--clear"]))]
            });
        }
        Ok(match mode.as_str() {
            "dust" => vec![("lcmask", if hard { s(&["--hard"]) } else { vec![] })],
            "soft" if hard => vec![("lcmask", s(&["--hard", "--no-dust"]))],
            "soft" => vec![],
            _ => vec![("lcmask", s(&["--clear"]))],
        })
    }

    fn masked(&mut self, src: &str, cmds: Cmds) -> String {
        if cmds.is_empty() {
            src.to_string()
        } else {
            self.chain(src.to_string(), cmds)
        }
    }

    fn search_args(&mut self, words: bool, hitchoice: bool, skip: &[&str]) -> Vec<String> {
        let mut out = Vec::new();
        let table: Vec<(&str, &str)> = CRITERIA
            .iter()
            .chain(if words { WORDS } else { &[] })
            .copied()
            .collect();
        for (v, r) in table {
            if !skip.contains(&v)
                && let Some(t) = self.a.take(v)
            {
                out.extend([r.to_string(), t]);
            }
        }
        for (v, r) in CRITERIA_FLAGS {
            if !skip.contains(v) && self.a.flag(v) {
                out.push(r.to_string());
            }
        }
        if hitchoice {
            if let Some(t) = self.a.take("maxhits") {
                out.extend(["--maxhits".to_string(), t]);
            }
            if self.a.flag("top_hits_only") {
                out.push("--top-hits-only".into());
            }
        }
        for v in IGNORED {
            if self.a.take(v).is_some() {
                self.caveat(format!("--{v} ignored, as vsearch ignores it"));
            }
        }
        out
    }

    fn hitfmt_scoring(&self) -> Vec<String> {
        let mut out = Vec::new();
        for (v, r) in HITFMT_SCORING {
            if let Some(t) = self.a.peek(v) {
                out.extend([r.to_string(), t.to_string()]);
            }
        }
        if self.a.has("n_mismatch") {
            out.push("--n-mismatch".into());
        }
        if self.cmd == "search_exact" {
            // every position of an exact hit scores a match
            out.push("--exact".into());
        }
        out
    }

    fn threads(&mut self) {
        let t = self.a.take("threads");
        if UNORDERED.contains(&self.cmd.as_str()) && t.as_deref() != Some("1") {
            self.caveat(
                "vsearch writes each query's records in the order its threads finish \
                 without --threads 1; rsearch writes them in input order"
                    .into(),
            );
        }
    }

    // the translation

    fn translate(mut self) -> R<Pipeline> {
        for v in ["quiet", "no_progress"] {
            self.a.take(v);
        }
        if !self.a.flag("notrunclabels") && !WHOLE_HEADERS.contains(&self.cmd.as_str()) {
            self.caveat(
                "rsearch keeps labels whole: vsearch cuts them at the first space or tab \
                 without --notrunclabels"
                    .into(),
            );
        }
        let gz = self.a.flag("gzip_decompress");
        let bz = self.a.flag("bzip2_decompress");
        if gz && bz {
            return err("--gzip_decompress with --bzip2_decompress");
        }
        if self.inp == "-" && (gz || bz) {
            return err("compressed standard input: give a file");
        }
        let inp = std::mem::take(&mut self.inp);
        self.inp = self.decompressed(inp);
        for v in ["reverse", "db"] {
            if let Some(path) = self.a.peek(v).map(str::to_string) {
                let d = self.decompressed(path);
                self.a.set(v, d);
            }
        }
        self.width = self.take_num("fasta_width", Some(80.0))?.unwrap();
        self.threads();
        let cmd = self.cmd.clone();
        if !HANDLED.contains(&cmd.as_str()) {
            return err(format!("--{cmd}: no rsearch translation"));
        }
        if self.inp.is_empty() {
            return err(format!("--{cmd} needs an input file"));
        }
        self.dispatch(&cmd)?;
        if self.a.has("log") && !self.a.used.contains("log") {
            self.a.take("log"); // a log of the run, not an output compared
        }
        let left = self.a.untaken();
        if !left.is_empty() {
            let names: Vec<String> = left.iter().map(|n| format!("--{n}")).collect();
            return err(format!(
                "{}: no rsearch translation for --{cmd}",
                names.join(", ")
            ));
        }
        Ok(self.p)
    }

    fn dispatch(&mut self, cmd: &str) -> R<()> {
        match cmd {
            "fastq_convert" => self.t_fastq_convert(),
            "fasta2fastq" => self.t_fasta2fastq(),
            "sff_convert" => self.t_sff_convert(),
            "udb2fasta" => self.t_udb2fasta(),
            "fastq_chars" => self.t_fastq_chars(),
            "fastq_stats" => self.t_fastq_stats(),
            "fastq_eestats" => self.t_fastq_eestats(),
            "fastq_eestats2" => self.t_fastq_eestats2(),
            "fastx_filter" | "fastq_filter" => self.t_fastx_filter(),
            "fastq_mergepairs" => self.t_fastq_mergepairs(),
            "fastq_join" => self.t_fastq_join(),
            "fastx_syncpairs" => self.t_fastx_syncpairs(),
            "fastx_uniques" => self.derep(&[], 1.0, "fastaout", Some("fastqout")),
            "derep_fulllength" => self.derep(&[], 32.0, "output", None),
            "derep_id" => self.derep(&["--same-label"], 32.0, "output", None),
            "derep_prefix" => self.derep(&["--prefix"], 32.0, "output", None),
            "derep_smallmem" => self.derep(&["--unsorted"], 1.0, "fastaout", None),
            "rereplicate" => self.t_rereplicate(),
            "fastx_revcomp" => self.t_fastx_revcomp(),
            "fastx_mask" => self.mask_command("fastaout", Some("fastqout")),
            "maskfasta" => self.mask_command("output", None),
            "cut" => self.t_cut(),
            "scramble" => self.t_scramble(),
            "fastx_getseq" | "fastx_getseqs" => self.getseq(false),
            "fastx_getsubseq" => self.getseq(true),
            "sortbylength" => self.sort("--length"),
            "sortbysize" => self.sort("--size"),
            "shuffle" => self.t_shuffle(),
            "fastx_subsample" => self.t_fastx_subsample(),
            "usearch_global" => self.t_usearch_global(),
            "search_exact" => self.t_search_exact(),
            "allpairs_global" => self.t_allpairs_global(),
            "orient" => self.t_orient(),
            "cluster_fast" => self.cluster(Some("--length"), false),
            "cluster_size" => self.cluster(Some("--size"), false),
            "cluster_smallmem" => self.cluster(None, false),
            "cluster_unoise" => self.cluster(Some("--size"), true),
            "uchime_denovo" => self.uchime(true, false, None),
            "uchime2_denovo" => self.uchime(true, true, None),
            "uchime3_denovo" => self.uchime(true, true, Some("16")),
            "uchime_ref" => self.uchime(false, false, None),
            "chimeras_denovo" => self.t_chimeras_denovo(),
            "sintax" => self.t_sintax(),
            "makeudb_usearch" => self.t_makeudb_usearch(),
            "udbinfo" => self.log_output(
                "udbinfo",
                vec!["-i".into(), self.inp.clone()],
                "udbinfo-log",
            ),
            "udbstats" => self.log_output(
                "udbstat",
                vec!["-i".into(), self.inp.clone()],
                "udbstats-log",
            ),
            _ => unreachable!(),
        }
    }

    /// A report that vsearch writes to --log: the tool's -o.
    fn log_output(&mut self, tool: &str, args: Vec<String>, comparator: &str) -> R<()> {
        if self.a.has("log") {
            let out = self.path("log");
            self.run(tool, args, None, Some(out.clone()), "out");
            self.output("log", &out, comparator);
        } else {
            self.run(tool, args, None, None, tool);
        }
        Ok(())
    }

    fn input_args(&self, rest: Vec<String>) -> Vec<String> {
        let mut args = vec!["-i".to_string(), self.inp.clone()];
        args.extend(rest);
        args
    }

    // format conversion

    fn t_fastq_convert(&mut self) -> R<()> {
        let mut args = self.qual_in();
        args.extend(self.qual_out());
        if self.a.flag("fastq_solexa") {
            args.push("--solexa".into());
        }
        let offset = self
            .a
            .peek("fastq_asciiout")
            .or(self.a.peek("fastq_ascii"))
            .unwrap_or("33")
            .to_string();
        let s = self.run("conv", self.input_args(args), None, None, "conv");
        let s = self.labels(&s, Some(&offset), &["sizein"]);
        self.fastq("fastqout", &s);
        self.a.take("sizein");
        Ok(())
    }

    fn t_fasta2fastq(&mut self) -> R<()> {
        let args = self.qual_out();
        let offset = self.opt_or("fastq_asciiout", "33");
        let s = self.run("conv", self.input_args(args), None, None, "conv");
        let s = self.labels(&s, Some(&offset), &[]);
        self.fastq("fastqout", &s);
        self.a.take("sizein");
        Ok(())
    }

    fn t_sff_convert(&mut self) -> R<()> {
        let mut args = if self.a.flag("sff_clip") {
            s(&["--clip"])
        } else {
            vec![]
        };
        args.extend(self.qual_out());
        let s = self.run("conv", self.input_args(args), None, None, "conv");
        let s = self.labels(&s, None, &[]);
        self.fastq("fastqout", &s);
        Ok(())
    }

    fn t_udb2fasta(&mut self) -> R<()> {
        if !self.label_args(None, &[]).is_empty() {
            return err("label options with --udb2fasta");
        }
        let inp = self.inp.clone();
        self.fasta("output", &inp, Some("udb"));
        Ok(())
    }

    // quality reports

    fn t_fastq_chars(&mut self) -> R<()> {
        let mut args = self.input_args(vec![]);
        if let Some(t) = self.a.take("fastq_tail") {
            args.extend(["--tail".to_string(), t]);
        }
        self.log_output("fqchar", args, "log-report")
    }

    fn t_fastq_stats(&mut self) -> R<()> {
        let q = self.qual_in();
        let args = self.input_args(q);
        self.log_output("fqstat", args, "log-report")
    }

    fn t_fastq_eestats(&mut self) -> R<()> {
        let q = self.qual_in();
        let s = self.run("eestat", self.input_args(q), None, None, "eestat");
        self.fastq("output", &s);
        Ok(())
    }

    fn t_fastq_eestats2(&mut self) -> R<()> {
        let mut args = self.input_args(vec![]);
        if let Some(t) = self.a.take("length_cutoffs") {
            let parts: Vec<&str> = t.split(',').collect();
            if parts.len() != 3 {
                return err("--length_cutoffs needs three values");
            }
            args.extend(s(&["--min-len", parts[0]]));
            if parts[1] != "*" {
                args.extend(s(&["--max-len", parts[1]]));
            }
            args.extend(s(&["--len-step", parts[2]]));
        }
        if let Some(t) = self.a.take("ee_cutoffs") {
            args.extend(["--max-ee".to_string(), t]);
        }
        args.extend(self.qual_in());
        let s = self.run("eegrid", args, None, None, "eegrid");
        self.fastq("output", &s);
        Ok(())
    }

    // trimming and filtering

    fn t_fastx_filter(&mut self) -> R<()> {
        let mut trim = self.mapped(&[
            ("fastq_stripleft", "--strip-left"),
            ("fastq_stripright", "--strip-right"),
            ("fastq_trunclen_keep", "--len"),
            ("fastq_truncqual", "--qual"),
            ("fastq_truncee", "--max-ee"),
            ("fastq_truncee_rate", "--max-ee-rate"),
        ]);
        let mut minlen = self.take_num("fastq_minlen", None)?;
        if self.a.has("fastq_trunclen") {
            if self.a.has("fastq_trunclen_keep") {
                return err("--fastq_trunclen with --fastq_trunclen_keep");
            }
            let n = self.take_num("fastq_trunclen", None)?.unwrap();
            trim.extend(["--len".to_string(), n.text.clone()]);
            minlen = Some(match minlen {
                Some(m) if m.v > n.v => m,
                _ => n,
            });
        }
        let mut filt = Vec::new();
        if let Some(m) = minlen {
            if m.v < 1.0 {
                return err(format!(
                    "--fastq_minlen {}: sqfilt keeps no empty read",
                    m.text
                ));
            }
            filt.extend(["--min-len".to_string(), m.text]);
        }
        filt.extend(self.mapped(&[
            ("fastq_maxlen", "--max-len"),
            ("fastq_maxns", "--max-ns"),
            ("fastq_maxee", "--max-ee"),
            ("fastq_maxee_rate", "--max-ee-rate"),
            ("fastq_minqual", "--min-qual"),
            ("minsize", "--min-size"),
            ("maxsize", "--max-size"),
        ]));
        let offset = self.opt_or("fastq_ascii", "33");
        let q = self.qual_in();
        let rev = self.a.take("reverse");
        let fmt = self.fmt(&self.inp.clone());
        let mut src = self.inp.clone();
        let mut rsrc = rev.clone();
        if !trim.is_empty() {
            let mut args = vec!["-i".to_string(), src];
            args.extend(trim.iter().cloned());
            args.extend(q.iter().cloned());
            src = self.run("sqtrim", args, None, None, "sqtrim");
            if let Some(r) = &rsrc {
                let mut args = vec!["-i".to_string(), r.clone()];
                args.extend(trim.iter().cloned());
                args.extend(q.iter().cloned());
                rsrc = Some(self.run("sqtrim", args, None, None, "sqtrim"));
            }
        }
        let kept = self.tmp("kept");
        let mut args = Vec::new();
        let (mut disc, mut rkept, mut rdisc) = (None, None, None);
        if let Some(r) = &rsrc {
            let k = self.tmp("rev_kept");
            args.extend([
                "--reverse".to_string(),
                r.clone(),
                "--rev-out".into(),
                k.clone(),
            ]);
            rkept = Some(k);
        }
        let wants = |a: &Args, o: &[&str]| o.iter().any(|x| a.has(x));
        if wants(&self.a, &["fastaout_discarded", "fastqout_discarded"]) {
            let d = self.tmp("discarded");
            args.extend(["--discarded".to_string(), d.clone()]);
            disc = Some(d);
        }
        if rev.is_some()
            && wants(
                &self.a,
                &["fastaout_discarded_rev", "fastqout_discarded_rev"],
            )
        {
            let d = self.tmp("rev_discarded");
            args.extend(["--rev-discarded".to_string(), d.clone()]);
            rdisc = Some(d);
        }
        args.extend(filt);
        args.extend(q);
        self.step("sqfilt", &src, args, Some(kept.clone()), "out");
        self.a.take("sizein");
        if disc.is_some() || rev.is_some() {
            self.no_labels("relabelling with several outputs is not translated", &[])?;
        }
        let kept = self.labels(&kept, Some(&offset), &[]);
        self.records(&kept, fmt, Some("fastaout"), Some("fastqout"));
        if let Some(d) = disc {
            self.records(
                &d,
                fmt,
                Some("fastaout_discarded"),
                Some("fastqout_discarded"),
            );
        }
        if let Some(rev) = rev {
            let rfmt = self.fmt(&rev);
            self.records(
                &rkept.unwrap(),
                rfmt,
                Some("fastaout_rev"),
                Some("fastqout_rev"),
            );
            if let Some(d) = rdisc {
                self.records(
                    &d,
                    rfmt,
                    Some("fastaout_discarded_rev"),
                    Some("fastqout_discarded_rev"),
                );
            }
        }
        Ok(())
    }

    // pairs

    fn t_fastq_mergepairs(&mut self) -> R<()> {
        let Some(mut rev) = self.a.take("reverse") else {
            return err("--fastq_mergepairs needs --reverse");
        };
        let q = self.qual_in();
        let offset = self.opt_or("fastq_ascii", "33");
        let mut fwd = self.inp.clone();
        let mut pre: Vec<String> = Vec::new();
        if let Some(t) = self.a.take("fastq_truncqual") {
            let trim = |file: &str| {
                let mut a = vec![
                    "-i".to_string(),
                    file.to_string(),
                    "--qual".into(),
                    t.clone(),
                ];
                a.extend(q.iter().cloned());
                a
            };
            let (af, ar) = (trim(&fwd), trim(&rev));
            fwd = self.run("sqtrim", af, None, None, "sqtrim");
            rev = self.run("sqtrim", ar, None, None, "sqtrim");
            pre.push("--fastq_truncqual".into());
        }
        let mut filt = Vec::new();
        for (v, r) in [
            ("fastq_minlen", "--min-len"),
            ("fastq_maxlen", "--max-len"),
            ("fastq_maxns", "--max-ns"),
        ] {
            if let Some(t) = self.a.take(v) {
                filt.extend([r.to_string(), t]);
                pre.push(format!("--{v}"));
            }
        }
        if !filt.is_empty() {
            let rkept = self.tmp("rev");
            let mut args = s(&["-i", &fwd, "--reverse", &rev, "--rev-out", &rkept]);
            args.extend(filt);
            args.extend(q.iter().cloned());
            fwd = self.run("sqfilt", args, None, None, "sqfilt");
            rev = rkept;
        }
        let mut margs = s(&["-i", &fwd, "--reverse", &rev]);
        margs.extend(self.mapped(&[
            ("fastq_minovlen", "--min-overlap"),
            ("fastq_maxdiffs", "--max-diffs"),
            ("fastq_maxdiffpct", "--max-diff-pct"),
            ("fastq_qminout", "--qmin-out"),
            ("fastq_qmaxout", "--qmax-out"),
        ]));
        if self.a.flag("fastq_allowmergestagger") {
            margs.push("--stagger".into());
        }
        self.a.take("fastq_nostagger");
        let post = self.mapped(&[
            ("fastq_maxee", "--max-ee"),
            ("fastq_minmergelen", "--min-len"),
            ("fastq_maxmergelen", "--max-len"),
        ]);
        let mut filters = pre.clone();
        filters.extend(post.iter().filter(|x| x.starts_with("--")).cloned());
        let mut side = Vec::new();
        for (opt, flag) in [
            ("notmerged_fwd", "--not-merged"),
            ("notmerged_rev", "--rev-not-merged"),
        ] {
            if self.a.has(&format!("fastaout_{opt}")) || self.a.has(&format!("fastqout_{opt}")) {
                if !filters.is_empty() {
                    return err(format!(
                        "--fast?out_{opt} with {}: the reads sqfilt discards are not written \
                         as not merged",
                        filters.join(", ")
                    ));
                }
                let f = self.tmp(opt);
                margs.extend([flag.to_string(), f.clone()]);
                side.push((opt, f));
            }
        }
        if self.a.has("eetabbedout") {
            if !post.is_empty() {
                return err("--eetabbedout with filters of merged reads");
            }
            margs.extend(["--ee-table".to_string(), self.path("eetabbedout")]);
        }
        let mut report = None;
        if self.a.has("log") {
            if !filters.is_empty() {
                self.a.take("log");
                self.caveat("the --log report is not compared: filters run outside mergpa".into());
            } else {
                let r = self.path("log");
                margs.extend(["--report".to_string(), r.clone()]);
                report = Some(r);
            }
        }
        margs.extend(q.iter().cloned());
        let mut merged = self.run("mergpa", margs, None, None, "merged");
        if self.a.has("eetabbedout") {
            let p = self.path("eetabbedout");
            self.output("eetabbedout", &p, "bytes");
        }
        if let Some(r) = report {
            self.output("log", &r, "log-report");
        }
        if !post.is_empty() {
            let mut args = post;
            args.extend(q);
            merged = self.run("sqfilt", args, Some(merged), None, "sqfilt");
        }
        self.a.take("sizein");
        let merged = self.labels(&merged, Some(&offset), &[]);
        self.records(&merged, Some("fastq"), Some("fastaout"), Some("fastqout"));
        for (opt, f) in side {
            self.records(
                &f,
                Some("fastq"),
                Some(&format!("fastaout_{opt}")),
                Some(&format!("fastqout_{opt}")),
            );
        }
        Ok(())
    }

    fn t_fastq_join(&mut self) -> R<()> {
        let Some(rev) = self.a.take("reverse") else {
            return err("--fastq_join needs --reverse");
        };
        let mut args = s(&["-i", &self.inp, "--reverse", &rev]);
        args.extend(self.mapped(&[
            ("join_padgap", "--pad"),
            ("join_padgapq", "--pad-qual"),
            ("fastq_ascii", "--ascii"),
        ]));
        let s = self.run("joinpa", args, None, None, "joinpa");
        self.a.take("sizein");
        let offset = self.opt_or("fastq_ascii", "33");
        let s = self.labels(&s, Some(&offset), &[]);
        self.records(&s, Some("fastq"), Some("fastaout"), Some("fastqout"));
        Ok(())
    }

    fn t_fastx_syncpairs(&mut self) -> R<()> {
        let Some(rev) = self.a.take("reverse") else {
            return err("--fastx_syncpairs needs --reverse");
        };
        let outs: Vec<String> = ["fwd", "rev", "orphans", "rev_orphans"]
            .iter()
            .map(|k| self.tmp(k))
            .collect();
        let mut args = s(&[
            "-i",
            &self.inp,
            "--reverse",
            &rev,
            "-o",
            &outs[0],
            "--rev-out",
            &outs[1],
            "--orphans",
            &outs[2],
            "--rev-orphans",
            &outs[3],
        ]);
        if let Some(t) = self.a.take("read_separators") {
            args.extend(["--separators".to_string(), t]);
        }
        self.run("syncpa", args, None, None, "syncpa");
        let (ffmt, rfmt) = (self.fmt(&self.inp.clone()), self.fmt(&rev));
        self.records(&outs[0], ffmt, Some("fastaout"), Some("fastqout"));
        self.records(&outs[1], rfmt, Some("fastaout_rev"), Some("fastqout_rev"));
        self.records(
            &outs[2],
            ffmt,
            Some("fastaout_orphans"),
            Some("fastqout_orphans"),
        );
        self.records(
            &outs[3],
            rfmt,
            Some("fastaout_orphans_rev"),
            Some("fastqout_orphans_rev"),
        );
        Ok(())
    }

    // dereplication

    fn derep(&mut self, mode: &[&str], min_default: f64, fa: &str, fq: Option<&str>) -> R<()> {
        let fmt = self.fmt(&self.inp.clone());
        let inp = self.inp.clone();
        let src = self.filtered(&inp, min_default, vec![])?;
        let mut args = s(mode);
        if let Some(st) = self.a.take("strand") {
            if st == "both" {
                args.push("--strand".into());
            } else if st != "plus" {
                return err(format!("--strand {st}"));
            }
        }
        let sizein = self.a.flag("sizein");
        let sizeout = self.a.flag("sizeout");
        if !sizein {
            args.push("--no-sizein".into());
            if src != self.inp {
                self.caveat(
                    "sqfilt checks ;size= annotations, which dereplication reads only with \
                     --sizein"
                        .into(),
                );
            }
        }
        if !sizeout {
            args.push("--no-sizeout".into());
        }
        for (v, r) in [("uc", "--uc"), ("tabbedout", "--tabbedout")] {
            if self.a.has(v) {
                args.extend([r.to_string(), self.path(v)]);
            }
        }
        args.extend(self.qual_in());
        args.extend(self.qual_out());
        if self.a.flag("fastq_qout_max") {
            args.push("--best-qual".into());
        }
        let mut st = self.step("derep", &src, args, None, "derep");
        for v in ["uc", "tabbedout"] {
            if self.a.has(v) {
                let p = self.path(v);
                self.output(v, &p, "bytes");
            }
        }
        let sizes = self.mapped(&[
            ("minuniquesize", "--min-size"),
            ("maxuniquesize", "--max-size"),
        ]);
        if !sizes.is_empty() {
            if !sizeout {
                return err("--minuniquesize or --maxuniquesize without --sizeout");
            }
            st = self.run("sqfilt", sizes, Some(st), None, "sqfilt");
        }
        if let Some(t) = self.a.take("topn") {
            st = self.run("topn", vec![t], Some(st), None, "topn");
        }
        let offset = self
            .a
            .peek("fastq_asciiout")
            .or(self.a.peek("fastq_ascii"))
            .unwrap_or("33")
            .to_string();
        let st = self.labels(&st, Some(&offset), &["sizeout"]);
        self.records(&st, fmt, Some(fa), fq);
        Ok(())
    }

    fn t_rereplicate(&mut self) -> R<()> {
        let extra = if self.a.flag("sizeout") {
            vec![]
        } else {
            s(&["--no-sizeout"])
        };
        self.a.take("sizein");
        let s = self.run("rerep", self.input_args(extra), None, None, "rerep");
        let s = self.labels(&s, None, &["sizeout"]);
        let fmt = self.fmt(&self.inp.clone());
        self.fasta("output", &s, fmt);
        Ok(())
    }

    // sequence manipulation

    fn t_fastx_revcomp(&mut self) -> R<()> {
        let s = self.run("revcomp", self.input_args(vec![]), None, None, "revcomp");
        self.a.take("sizein");
        let offset = self.opt_or("fastq_ascii", "33");
        let s = self.labels(&s, Some(&offset), &[]);
        let fmt = self.fmt(&self.inp.clone());
        self.records(&s, fmt, Some("fastaout"), Some("fastqout"));
        Ok(())
    }

    fn mask_command(&mut self, fa: &str, fq: Option<&str>) -> R<()> {
        let mode = self.take_or("qmask", "dust");
        let hard = self.a.flag("hardmask");
        let cmds: Cmds = match mode.as_str() {
            "dust" => vec![("lcmask", if hard { s(&["--hard"]) } else { vec![] })],
            "soft" if hard => vec![("lcmask", s(&["--hard", "--no-dust"]))],
            "soft" | "none" => vec![],
            _ => return err(format!("--qmask {mode}")),
        };
        let mut pct = Vec::new();
        for (v, which) in [("min_unmasked_pct", "min"), ("max_unmasked_pct", "max")] {
            if self.a.has(v) {
                if mode == "none" {
                    return err(format!("--{v} with --qmask none"));
                }
                let t = self.a.take(v).unwrap();
                let kind = if hard { "non-n" } else { "upper" };
                pct.extend([format!("--{which}-{kind}-pct"), t]);
            }
        }
        let inp = self.inp.clone();
        let has_cmds = !cmds.is_empty();
        let mut st = self.masked(&inp, cmds);
        if let Some(args) = self.length_filter(1.0, pct, "input")? {
            st = if has_cmds {
                self.run("sqfilt", args, Some(st), None, "sqfilt")
            } else {
                let mut all = vec!["-i".to_string(), st];
                all.extend(args);
                self.run("sqfilt", all, None, None, "sqfilt")
            };
        }
        self.a.take("sizein");
        let offset = self.opt_or("fastq_ascii", "33");
        let st = self.labels(&st, Some(&offset), &[]);
        let fmt = self.fmt(&inp);
        self.records(&st, fmt, Some(fa), fq);
        Ok(())
    }

    fn t_cut(&mut self) -> R<()> {
        if !self.a.has("cut_pattern") {
            return err("--cut needs --cut_pattern");
        }
        self.no_labels(
            "relabelling the four outputs is not translated",
            &["sizein"],
        )?;
        self.a.take("sizein");
        let cut = self.tmp("cut");
        let pattern = self.a.take("cut_pattern").unwrap();
        let mut args = s(&["-i", &self.inp, "--pattern", &pattern]);
        let mut outs = Vec::new();
        // the other outputs only where vsearch writes them
        for (key, opt, flag) in [
            ("rev", "fastaout_rev", "--rev-out"),
            ("uncut", "fastaout_discarded", "--uncut"),
            ("rev_uncut", "fastaout_discarded_rev", "--rev-uncut"),
        ] {
            if self.a.has(opt) {
                let f = self.tmp(key);
                args.extend([flag.to_string(), f.clone()]);
                outs.push((opt, f));
            }
        }
        self.run("redig", args, None, Some(cut.clone()), "out");
        let fmt = self.fmt(&self.inp.clone());
        self.fasta("fastaout", &cut, fmt);
        for (opt, f) in outs {
            self.fasta(opt, &f, fmt);
        }
        Ok(())
    }

    fn seed(&mut self) -> Vec<String> {
        match self.a.take("randseed") {
            Some(t) if t != "0" => vec![format!("--seed={t}")],
            _ => {
                self.caveat(
                    "no --randseed: vsearch and rsearch draw different random numbers".into(),
                );
                vec![]
            }
        }
    }

    fn t_scramble(&mut self) -> R<()> {
        let mut args = self.seed();
        if let Some(t) = self.a.take("scramble_kmer") {
            args.extend(["--kmer".to_string(), t]);
        }
        let s = self.run("scrmbl", self.input_args(args), None, None, "scrmbl");
        self.a.take("sizein");
        let offset = self.opt_or("fastq_ascii", "33");
        let s = self.labels(&s, Some(&offset), &[]);
        let fmt = self.fmt(&self.inp.clone());
        self.records(&s, fmt, Some("fastaout"), Some("fastqout"));
        Ok(())
    }

    // retrieval

    fn select_args(&mut self) -> Vec<String> {
        let mut args = self.mapped(&[
            ("label", "--label"),
            ("labels", "--labels"),
            ("label_word", "--word"),
            ("label_words", "--words"),
            ("label_field", "--field"),
        ]);
        if self.a.flag("label_substr_match") {
            args.push("--substr".into());
        }
        args
    }

    fn getseq(&mut self, region: bool) -> R<()> {
        let sel = self.select_args();
        if region && !sel.iter().any(|a| a == "--label") {
            return err("--fastx_getsubseq needs --label");
        }
        let kept = self.tmp("kept");
        let mut args = self.input_args(sel);
        let mut disc = None;
        if self.a.has("notmatched") || self.a.has("notmatchedfq") {
            let d = self.tmp("notmatched");
            args.extend(["--discarded".to_string(), d.clone()]);
            disc = Some(d);
            self.no_labels(
                "relabelling with --notmatched is not translated",
                &["sizein"],
            )?;
        }
        self.run("getseq", args, None, Some(kept.clone()), "out");
        let mut st = kept;
        if region {
            let sub = self.mapped(&[("subseq_start", "--start"), ("subseq_end", "--end")]);
            st = self.run("subseq", sub, Some(st), None, "subseq");
        }
        self.a.take("sizein");
        let fmt = self.fmt(&self.inp.clone());
        let offset = self.opt_or("fastq_ascii", "33");
        let st = self.labels(&st, Some(&offset), &[]);
        self.records(&st, fmt, Some("fastaout"), Some("fastqout"));
        if let Some(d) = disc {
            self.fasta("notmatched", &d, fmt);
            self.fastq("notmatchedfq", &d);
        }
        Ok(())
    }

    // sorting and sampling

    fn sort(&mut self, key: &str) -> R<()> {
        let extra = if key == "--size" {
            self.mapped(&[("minsize", "--min-size"), ("maxsize", "--max-size")])
        } else {
            vec![]
        };
        let inp = self.inp.clone();
        let st = self.filtered(&inp, 1.0, extra)?;
        let mut st = self.chain(st, vec![("sortby", s(&[key]))]);
        if let Some(t) = self.a.take("topn") {
            st = self.run("topn", vec![t], Some(st), None, "topn");
        }
        self.a.take("sizein");
        let st = self.labels(&st, None, &[]);
        let fmt = self.fmt(&inp);
        self.fasta("output", &st, fmt);
        Ok(())
    }

    fn t_shuffle(&mut self) -> R<()> {
        let inp = self.inp.clone();
        let st = self.filtered(&inp, 1.0, vec![])?;
        let seed = self.seed();
        let mut st = self.step("shufl", &st, seed, None, "shufl");
        if let Some(t) = self.a.take("topn") {
            st = self.run("topn", vec![t], Some(st), None, "topn");
        }
        self.a.take("sizein");
        let st = self.labels(&st, None, &[]);
        let fmt = self.fmt(&inp);
        self.fasta("output", &st, fmt);
        Ok(())
    }

    fn t_fastx_subsample(&mut self) -> R<()> {
        let inp = self.inp.clone();
        let src = self.filtered(&inp, 1.0, vec![])?;
        let mut args = self.mapped(&[("sample_size", "--count"), ("sample_pct", "--pct")]);
        let sizein = self.a.flag("sizein");
        let sizeout = self.a.flag("sizeout");
        if sizein != sizeout {
            return err(
                "--sizein without --sizeout, or the other way round: subsmp --by-size does both",
            );
        }
        if sizein {
            args.push("--by-size".into());
        }
        if self.a.flag("allow_fewer") {
            args.push("--allow-fewer".into());
        }
        args.extend(self.seed());
        let mut disc = None;
        if self.a.has("fastaout_discarded") || self.a.has("fastqout_discarded") {
            let d = self.tmp("discarded");
            args.extend(["--discarded".to_string(), d.clone()]);
            disc = Some(d);
        }
        self.no_labels("relabelling subsamples is not translated", &[])?;
        let kept = self.step("subsmp", &src, args, None, "subsmp");
        let fmt = self.fmt(&inp);
        self.records(&kept, fmt, Some("fastaout"), Some("fastqout"));
        if let Some(d) = disc {
            self.records(
                &d,
                fmt,
                Some("fastaout_discarded"),
                Some("fastqout_discarded"),
            );
        }
        Ok(())
    }

    // searches: hits and what is made from them

    /// hitfmt, getseq, otutab and lcatax on the hits. `hits.get(kind)` gives
    /// the hit file: 'all' every hit the search keeps, 'top' its first hit
    /// a query (vsearch's --uc without --uc_allhits), 'every' without
    /// --maxhits or --top_hits_only (for --dbmatched, which counts every
    /// accepted hit).
    fn hit_outputs(&mut self, hits: &mut Hits, query: &str, db: &str) -> R<()> {
        let scoring = self.hitfmt_scoring();
        let nohits = if self.a.flag("output_no_hits") {
            s(&["--output-no-hits"])
        } else {
            vec![]
        };
        if self.a.has("uc") {
            let kind = if self.a.flag("uc_allhits") {
                "all"
            } else {
                "top"
            };
            let h = hits.get(self, kind);
            // vsearch's std::sort leaves hits that tie on its sort key (the
            // two strands of one target) in either order, where rsearch
            // writes the plus strand first. The comparators accept both,
            // reading the queries to know which are their own reverse
            // complement
            let ties = format!("uc-ties:{}", self.inp);
            self.output("uc", &h, &ties);
        }
        self.a.take("uc_allhits");
        let fmts: [(&str, &[&str], &str); 7] = [
            ("alnout", &["--alignments"], "alnout"),
            ("blast6out", &["--blast6"], "blast6-ties"),
            ("userout", &[], "bytes"),
            ("samout", &["--sam"], "sam-ties"),
            ("fastapairs", &["--fasta-pairs"], "bytes"),
            ("qsegout", &["--qsegs"], "bytes"),
            ("tsegout", &["--tsegs"], "bytes"),
        ];
        for (opt, base, comp) in fmts {
            if !self.a.has(opt) {
                continue;
            }
            let mut args = s(base);
            if opt == "userout" {
                let Some(fields) = self.a.take("userfields") else {
                    return err("--userout without --userfields");
                };
                args = vec!["--userout".into(), fields];
            }
            if opt == "alnout"
                && let Some(t) = self.a.take("rowlen")
            {
                args.extend(["--rowlen".to_string(), t]);
            }
            if opt == "samout" && self.a.flag("samheader") {
                args.push("--sam-header".into());
            }
            if opt == "fastapairs" {
                self.unwrapped_only(opt)?;
            }
            let h = hits.get(self, "all");
            let mut all = s(&[&h, "-i", query, "--db", db]);
            all.extend(args);
            all.extend(scoring.iter().cloned());
            all.extend(nohits.iter().cloned());
            let out = self.run("hitfmt", all, None, None, opt);
            if opt == "qsegout" || opt == "tsegout" {
                self.fasta(opt, &out, Some("fasta"));
            } else if comp.ends_with("-ties") {
                let comp = format!("{comp}:{}", self.inp);
                self.output(opt, &out, &comp);
            } else {
                self.output(opt, &out, comp);
            }
        }
        for (opt, kind, extra, src) in [
            ("matched", "all", &[][..], query),
            ("notmatched", "all", &["--invert"][..], query),
            ("dbmatched", "every", &["--targets"][..], db),
            ("dbnotmatched", "every", &["--targets", "--invert"][..], db),
        ] {
            if !self.a.has(opt) {
                continue;
            }
            if opt.starts_with("db") && self.a.has("weak_id") {
                self.caveat(
                    "weak hits count for getseq --targets, which vsearch does not count for \
                     --dbmatched"
                        .into(),
                );
            }
            let h = hits.get(self, kind);
            let mut args = s(&["--hits", &h, "-i", src]);
            args.extend(s(extra));
            let st = self.run("getseq", args, None, None, opt);
            let fmt = if src != db {
                self.fmt(src)
            } else {
                Some("fasta")
            };
            self.fasta(opt, &st, fmt);
        }
        for (opt, extra, comp) in [
            ("otutabout", &[][..], "bytes"),
            ("mothur_shared_out", &["--mothur"][..], "bytes"),
            ("biomout", &["--biom"][..], "biom"),
        ] {
            if self.a.has(opt) {
                let h = hits.get(self, "all");
                let mut args = s(&[&h, "--db", db]);
                args.extend(s(extra));
                let out = self.run("otutab", args, None, None, opt);
                self.output(opt, &out, comp);
            }
        }
        if self.a.has("lcaout") {
            let h = hits.get(self, "all");
            let mut args = s(&[&h, "--db", db]);
            if let Some(t) = self.a.take("lca_cutoff") {
                args.extend(["--cutoff".to_string(), t]);
            }
            let out = self.run("lcatax", args, None, None, "lcaout");
            self.output("lcaout", &out, "bytes");
        }
        Ok(())
    }

    fn no_search_labels(&mut self) -> R<()> {
        self.a.take("sizein");
        self.no_labels(
            &format!("--{} relabels nothing rsearch writes", self.cmd),
            &[],
        )
    }

    /// The database as searched and as written (hitfmt, getseq): sqfilt
    /// for the length filters, lcmask for --dbmask. A UDB file is taken as
    /// it is.
    fn database(&mut self, min_default: f64, mask_style: Option<&str>) -> R<(String, String)> {
        let Some(db) = self.a.take("db") else {
            return err(format!("--{} needs --db", self.cmd));
        };
        if self.fmt(&db) == Some("udb") {
            self.lengths(min_default)?;
            if mask_style.is_some() {
                self.a.take("dbmask");
                self.a.take("hardmask");
            }
            return Ok((db.clone(), db));
        }
        let (mn, mx) = self.lengths(min_default)?;
        let filt = if mn.v <= 0.0 {
            None
        } else {
            let mut f = if mn.v == 1.0 {
                vec![]
            } else {
                vec!["--min-len".into(), mn.text.clone()]
            };
            f.extend(["--max-len".to_string(), mx.text.clone()]);
            Some(f)
        };
        if mn.v <= 0.0 {
            self.no_max_caveat(&mn, &mx, "database");
        }
        let mode = mask_style.map(|_| self.opt_or("dbmask", "dust"));
        let cmds = match mask_style {
            Some(style) => self.masking("dbmask", style)?,
            None => vec![],
        };
        let build = |t: &mut Translator, cmds: Cmds| {
            let mut st = t.masked(&db, cmds);
            if let Some(f) = &filt {
                st = if st == db {
                    let mut all = vec!["-i".to_string(), st];
                    all.extend(f.iter().cloned());
                    t.run("sqfilt", all, None, None, "db")
                } else {
                    t.run("sqfilt", f.clone(), Some(st), None, "db")
                };
            }
            st
        };
        let searched = build(self, cmds);
        let shown = if mode.as_deref() != Some("none") {
            searched.clone()
        } else {
            build(self, vec![])
        };
        Ok((searched, shown))
    }

    fn queries(&mut self) -> R<(String, String)> {
        let mode = self.opt_or("qmask", "dust");
        let cmds = self.masking("qmask", "standard")?;
        let inp = self.inp.clone();
        let searched = self.masked(&inp, cmds);
        let shown = if mode == "none" {
            inp
        } else {
            searched.clone()
        };
        Ok((searched, shown))
    }

    fn t_usearch_global(&mut self) -> R<()> {
        if !self.a.has("id") {
            return err("--usearch_global needs --id");
        }
        if self.opt_or("qmask", "dust") == "dust" && self.a.peek("strand") == Some("both") {
            self.caveat(
                "--qmask dust with --strand both: vsearch dusts each strand of a query on its \
                 own, where lcmask masks the plus strand and dbsrch reverse-complements it"
                    .into(),
            );
        }
        let (q, qshown) = self.queries()?;
        let (db, dbshown) = self.database(32.0, Some("standard"))?;
        let mut choice = Vec::new();
        if let Some(t) = self.a.take("maxhits") {
            choice.extend(["--maxhits".to_string(), t]);
        }
        if self.a.flag("top_hits_only") {
            choice.push("--top-hits-only".into());
        }
        let mut base = s(&["--db", &db]);
        base.extend(self.search_args(true, false, &[]));
        let mut hits = Hits::new("dbsrch", &q, base, choice);
        self.no_search_labels()?;
        self.hit_outputs(&mut hits, &qshown, &dbshown)
    }

    fn t_search_exact(&mut self) -> R<()> {
        if self.a.flag("hardmask") {
            return err("--hardmask with --search_exact");
        }
        let mode = self.take_or("qmask", "dust");
        self.a.take("dbmask");
        let (db, _) = self.database(1.0, None)?;
        let mut dbshown = db.clone();
        let mut qshown = self.inp.clone();
        let shows = [
            "alnout",
            "blast6out",
            "userout",
            "samout",
            "fastapairs",
            "qsegout",
            "tsegout",
            "matched",
            "notmatched",
            "dbmatched",
            "dbnotmatched",
        ]
        .iter()
        .any(|o| self.a.has(o));
        let dbmode = self.opt_or("dbmask", "dust");
        if shows && (mode == "dust" || dbmode == "dust") {
            self.caveat(
                "the sequences --search_exact writes as vsearch masked them: unverified".into(),
            );
            if mode == "dust" {
                qshown = self.chain(self.inp.clone(), vec![("lcmask", vec![])]);
            }
            if dbmode == "dust" {
                dbshown = self.chain(db.clone(), vec![("lcmask", vec![])]);
            }
        }
        let mut base = s(&["--db", &db]);
        if let Some(t) = self.a.take("strand") {
            base.extend(["--strand".to_string(), t]);
        }
        if self.a.flag("self") {
            base.push("--self".into());
        }
        base.extend(self.mapped(&[
            ("maxqsize", "--max-qsize"),
            ("mintsize", "--min-tsize"),
            ("minsizeratio", "--min-size-ratio"),
            ("maxsizeratio", "--max-size-ratio"),
        ]));
        let choice = match self.a.peek("maxhits") {
            Some(t) if t != "0" => s(&["--maxhits", t]),
            _ => vec![],
        };
        self.a.take("maxhits");
        self.a.take("top_hits_only"); // every exact hit is a top hit
        let inp = self.inp.clone();
        let mut hits = Hits::new("exsrch", &inp, base, choice);
        self.no_search_labels()?;
        self.hit_outputs(&mut hits, &qshown, &dbshown)?;
        // the scores reach only hitfmt (raw), which hit_outputs gave them
        self.a.take("match");
        self.a.take("mismatch");
        Ok(())
    }

    fn t_allpairs_global(&mut self) -> R<()> {
        let inp = self.inp.clone();
        let src = self.filtered(&inp, 1.0, vec![])?;
        let mode = self.opt_or("qmask", "dust");
        let cmds = self.masking("qmask", "standard")?;
        let f = self.masked(&src, cmds);
        let shown = if mode == "none" { src } else { f.clone() };
        let mut args = Vec::new();
        if self.a.flag("acceptall") {
            args.push("--accept-all".into());
        }
        for v in ["maxaccepts", "maxrejects", "wordlength", "minwordmatches"] {
            self.a.take(v); // vsearch's do nothing: allpr aligns every pair
        }
        args.extend(self.search_args(false, true, &[]));
        let h = self.step("allpr", &f, args, None, "uc");
        self.no_search_labels()?;
        let mut hits = Hits::fixed(h);
        self.hit_outputs(&mut hits, &shown, &shown)
    }

    fn t_orient(&mut self) -> R<()> {
        let qmode = self.take_or("qmask", "dust");
        if !["dust", "soft", "none"].contains(&qmode.as_str()) {
            return err(format!("--qmask {qmode}"));
        }
        let mut q = self.inp.clone();
        if qmode == "none" {
            for o in ["fastaout", "fastqout", "notmatched"] {
                if self.a.has(o) {
                    return err(format!(
                        "--{o} with --qmask none: lcmask --clear upper-cases the queries vsearch \
                         writes as they are"
                    ));
                }
            }
            q = self.chain(q, vec![("lcmask", s(&["--clear"]))]);
        }
        let Some(mut db) = self.a.take("db") else {
            return err("--orient needs --db");
        };
        if self.fmt(&db) == Some("udb") {
            self.lengths(1.0)?;
            self.a.take("dbmask");
            self.a.take("hardmask");
        } else {
            let mode = self.take_or("dbmask", "dust");
            let hard = self.a.flag("hardmask");
            let cmds: Cmds = match (mode.as_str(), hard) {
                ("dust", false) => vec![("lcmask", vec![])],
                ("dust", true) => vec![("lcmask", s(&["--hard"]))],
                ("soft", false) => vec![],
                ("soft", true) => vec![("lcmask", s(&["--hard", "--no-dust"]))],
                ("none", _) => vec![("lcmask", s(&["--clear"]))],
                _ => return err(format!("--dbmask {mode}")),
            };
            let mut st = db.clone();
            if let Some(args) = self.length_filter(1.0, vec![], "database")? {
                let mut all = vec!["-i".to_string(), db.clone()];
                all.extend(args);
                st = self.run("sqfilt", all, None, None, "db");
            }
            db = if cmds.is_empty() {
                st
            } else {
                self.chain(st, cmds)
            };
        }
        let mut args = s(&["--db", &db]);
        if let Some(t) = self.a.take("wordlength") {
            args.extend(["--word-len".to_string(), t]);
        }
        let mut notm = None;
        if self.a.has("notmatched") {
            let n = self.tmp("notmatched");
            args.extend(["--not-oriented".to_string(), n.clone()]);
            notm = Some(n);
        }
        if self.a.has("tabbedout") {
            args.extend(["--tabbedout".to_string(), self.path("tabbedout")]);
        }
        if self.a.has("log") {
            args.extend(["--report".to_string(), self.path("log")]);
        }
        let out = self.step("dbornt", &q, args, None, "oriented");
        if self.a.has("tabbedout") {
            let p = self.path("tabbedout");
            self.output("tabbedout", &p, "bytes");
        }
        if self.a.has("log") {
            let p = self.path("log");
            self.output("log", &p, "orient-log");
        }
        self.no_search_labels()?;
        let fmt = self.fmt(&self.inp.clone());
        self.records(&out, fmt, Some("fastaout"), Some("fastqout"));
        if self.a.has("notmatched") {
            let n = notm.unwrap();
            if fmt == Some("fastq") {
                self.output("notmatched", &n, "bytes");
            } else {
                self.fasta("notmatched", &n, fmt);
            }
        }
        Ok(())
    }

    // clustering

    fn cluster(&mut self, sort: Option<&str>, unoise: bool) -> R<()> {
        let mut extra = Vec::new();
        if unoise {
            let m = self.take_or("minsize", "8");
            if m != "1" {
                extra = vec!["--min-size".to_string(), m];
            }
        }
        let has_extra = !extra.is_empty();
        let inp = self.inp.clone();
        let src = self.filtered(&inp, 32.0, extra)?;
        if has_extra && src == inp {
            return err("--minsize with --minseqlength 0");
        }
        let mut cmds = self.masking("qmask", "standard")?;
        if let Some(key) = sort {
            cmds.push(("sortby", s(&[key])));
        }
        let st = self.masked(&src, cmds);
        let mut args = Vec::new();
        if unoise {
            self.a.take("id"); // vsearch ignores --id for --cluster_unoise
            args.extend(self.mapped(&[("weak_id", "--id"), ("unoise_alpha", "--alpha")]));
        } else if !self.a.has("id") {
            return err(format!("--{} needs --id", self.cmd));
        }
        if self.cmd == "cluster_fast" && !self.a.has("maxrejects") {
            args.extend(s(&["--maxrejects", "8"]));
        }
        if self.cmd == "cluster_smallmem" && !self.a.flag("usersort") {
            self.caveat("without --usersort vsearch requires the input sorted by length".into());
        }
        let skip: &[&str] = if unoise { &["weak_id", "id"] } else { &[] };
        args.extend(self.search_args(true, false, skip));
        if self.a.flag("sizeorder") {
            args.push("--size-order".into());
        }
        let sizein = self.a.flag("sizein");
        let sizeout = self.a.flag("sizeout");
        let mut sz = Vec::new();
        if !sizein {
            sz.push("--no-sizein".to_string());
        }
        if !sizeout {
            sz.push("--no-sizeout".to_string());
        }
        args.extend(sz.iter().cloned());
        let mut hits = None;
        if ["uc", "clusters", "msaout", "consout", "profile"]
            .iter()
            .any(|o| self.a.has(o))
        {
            let h = self.tmp("uc");
            args.extend(["--hits".to_string(), h.clone()]);
            hits = Some(h);
        }
        let tool = if unoise { "unoise" } else { "clust" };
        let cent = self.step(tool, &st, args, None, "centroids");
        self.no_labels(
            &format!(
                "relabelling the centroids of --{} is not translated",
                self.cmd
            ),
            &["sizeout"],
        )?;
        for o in [
            "alnout",
            "blast6out",
            "userout",
            "samout",
            "fastapairs",
            "qsegout",
            "tsegout",
            "matched",
            "notmatched",
            "otutabout",
            "mothur_shared_out",
            "biomout",
            "centroid_sizeout",
            "clusterout_id",
            "clusterout_sort",
            "cons_truncate",
        ] {
            if self.a.has(o) {
                return err(format!("--{o} with --{}: no rsearch translation", self.cmd));
            }
        }
        let fmt = self.fmt(&inp);
        self.fasta("centroids", &cent, fmt);
        let hits = hits.unwrap_or_default();
        if self.a.has("uc") {
            let comp = if self.cmd == "cluster_fast" {
                "uc-fast"
            } else {
                "bytes"
            };
            self.output("uc", &hits, comp);
        }
        if self.a.has("clusters") {
            if self.width.v != 0.0 && fmt != Some("fastq") {
                return err(format!(
                    "--clusters with --fasta_width {}: clspl writes records unwrapped (give \
                     --fasta_width 0)",
                    self.width.text
                ));
            }
            if sizeout {
                self.caveat(
                    "--clusters with --sizeout: vsearch rewrites the abundances clspl writes as \
                     they are"
                        .into(),
                );
            }
            let prefix = self.path("clusters");
            self.run(
                "clspl",
                s(&["--hits", &hits, "-i", &st, "--prefix", &prefix]),
                None,
                None,
                "clspl",
            );
            self.output("clusters", &prefix, "cluster-files");
        }
        for (opt, flag) in [
            ("msaout", &[][..]),
            ("consout", &["--consensus"][..]),
            ("profile", &["--profile"][..]),
        ] {
            if !self.a.has(opt) {
                continue;
            }
            if opt != "consout" {
                self.unwrapped_only(opt)?;
            }
            let mut args = s(&["--hits", &hits, "-i", &st]);
            args.extend(s(flag));
            args.extend(sz.iter().cloned());
            let out = self.run("clmsa", args, None, None, opt);
            if opt == "consout" {
                self.fasta(opt, &out, Some("fasta"));
            } else {
                self.output(opt, &out, "bytes");
            }
        }
        Ok(())
    }

    // chimeras

    fn uchime(&mut self, denovo: bool, uchime2: bool, abskew: Option<&str>) -> R<()> {
        let (st, mut args) = if denovo {
            let inp = self.inp.clone();
            let src = self.filtered(&inp, 1.0, vec![])?;
            let mut cmds = self.masking("qmask", "standard")?;
            cmds.push(("sortby", s(&["--size"])));
            (self.masked(&src, cmds), vec![])
        } else {
            // vsearch masks no query of --uchime_ref: dust and soft both skip
            // the words of its lower case as read, which uchime does, and it
            // writes the queries as read
            let mode = self.take_or("qmask", "dust");
            let inp = self.inp.clone();
            let q = match mode.as_str() {
                "dust" | "soft" => inp,
                "none" => {
                    self.caveat(
                        "--qmask none: lcmask --clear upper-cases the queries, which --uchime_ref \
                         writes as read"
                            .into(),
                    );
                    self.chain(inp, vec![("lcmask", s(&["--clear"]))])
                }
                m => return err(format!("--qmask {m}: not dust, soft or none")),
            };
            let (db, _) = self.database(1.0, Some("standard"))?;
            self.a.take("abskew"); // used de novo only
            if self.a.peek("strand") == Some("plus") {
                self.a.take("strand");
            }
            (q, s(&["--db", &db]))
        };
        if uchime2 {
            args.push("--uchime2".into());
            for v in ["minh", "mindiv", "mindiffs"] {
                self.a.take(v); // vsearch ignores them for UCHIME2
            }
        }
        if denovo {
            let skew = self.a.take("abskew").or(abskew.map(str::to_string));
            if let Some(k) = skew {
                args.extend(["--abskew".to_string(), k]);
            }
        }
        for v in ["minh", "mindiv", "mindiffs", "xn", "dn"] {
            if let Some(t) = self.a.take(v) {
                args.extend([format!("--{v}"), t]);
            }
        }
        self.a.take("sizein");
        self.no_labels(
            "uchime writes its records unchanged, relabelling is relab's",
            &[],
        )?;
        let mut outs = Vec::new();
        for (v, r) in [("chimeras", "--chimeras"), ("borderline", "--borderline")] {
            if self.a.has(v) {
                let f = self.tmp(v);
                args.extend([r.to_string(), f.clone()]);
                outs.push((v, f));
            }
        }
        for (v, r) in [("uchimeout", "--table"), ("uchimealns", "--alignments")] {
            if self.a.has(v) {
                args.extend([r.to_string(), self.path(v)]);
            }
        }
        if let Some(t) = self.a.take("alignwidth") {
            args.extend(["--width".to_string(), t]);
        }
        let non = self.step("uchime", &st, args, None, "nonchimeras");
        for v in ["uchimeout", "uchimealns"] {
            if self.a.has(v) {
                let p = self.path(v);
                self.output(v, &p, "bytes");
            }
        }
        let fmt = self.fmt(&self.inp.clone());
        self.fasta("nonchimeras", &non, fmt);
        for (v, f) in outs {
            self.fasta(v, &f, fmt);
        }
        Ok(())
    }

    fn t_chimeras_denovo(&mut self) -> R<()> {
        let inp = self.inp.clone();
        let src = self.filtered(&inp, 1.0, vec![])?;
        let mut cmds = self.masking("qmask", "standard")?;
        cmds.push(("sortby", s(&["--size"])));
        let st = self.masked(&src, cmds);
        let mut args = self.mapped(&[
            ("chimeras_parts", "--parts"),
            ("chimeras_parents_max", "--max-parents"),
            ("chimeras_length_min", "--min-segment"),
            ("chimeras_diff_pct", "--diff-pct"),
            ("abskew", "--abskew"),
            ("alignwidth", "--width"),
        ]);
        self.a.take("sizein");
        self.no_labels(
            "mchime writes its records unchanged, relabelling is relab's",
            &[],
        )?;
        let mut chim = None;
        if self.a.has("chimeras") {
            let c = self.tmp("chimeras");
            args.extend(["--chimeras".to_string(), c.clone()]);
            chim = Some(c);
        }
        for (v, r) in [("tabbedout", "--table"), ("alnout", "--alignments")] {
            if self.a.has(v) {
                args.extend([r.to_string(), self.path(v)]);
            }
        }
        let non = self.step("mchime", &st, args, None, "nonchimeras");
        for v in ["tabbedout", "alnout"] {
            if self.a.has(v) {
                let p = self.path(v);
                self.output(v, &p, "bytes");
            }
        }
        let fmt = self.fmt(&inp);
        self.fasta("nonchimeras", &non, fmt);
        if let Some(c) = chim {
            self.fasta("chimeras", &c, fmt);
        }
        Ok(())
    }

    // taxonomy

    fn t_sintax(&mut self) -> R<()> {
        let Some(mut db) = self.a.take("db") else {
            return err("--sintax needs --db");
        };
        if self.fmt(&db) != Some("udb") {
            let mode = self.take_or("dbmask", "dust");
            if !["dust", "soft", "none"].contains(&mode.as_str()) {
                return err(format!("--dbmask {mode}"));
            }
            let mut st = db.clone();
            if let Some(args) = self.length_filter(32.0, vec![], "database")? {
                let mut all = vec!["-i".to_string(), db.clone()];
                all.extend(args);
                st = self.run("sqfilt", all, None, None, "db");
            }
            if mode == "none" {
                st = self.chain(st, vec![("lcmask", s(&["--clear"]))]);
            }
            db = st;
        } else {
            self.lengths(32.0)?;
            self.a.take("dbmask");
        }
        let q = self.chain(self.inp.clone(), vec![("lcmask", s(&["--clear"]))]);
        let mut args = s(&["--db", &db]);
        args.extend(self.mapped(&[
            ("wordlength", "--word-len"),
            ("strand", "--strand"),
            ("sintax_cutoff", "--cutoff"),
        ]));
        if self.a.flag("sintax_random") {
            args.push("--random".into());
        }
        args.extend(self.seed());
        if self.a.has("log") {
            args.extend(["--report".to_string(), self.path("log")]);
        }
        let out = self.run("sintax", args, Some(q), None, "sintax");
        self.fastq("tabbedout", &out);
        if self.a.has("log") {
            let p = self.path("log");
            self.output("log", &p, "sintax-log");
        }
        Ok(())
    }

    // UDB files

    fn t_makeudb_usearch(&mut self) -> R<()> {
        let inp = self.inp.clone();
        let src = self.filtered(&inp, 32.0, vec![])?;
        let cmds = self.masking("dbmask", "udb")?;
        let st = self.masked(&src, cmds);
        let args = self.mapped(&[("wordlength", "--word-len")]);
        let out = self.step("mkudb", &st, args, None, "udb");
        self.fastq("output", &out);
        Ok(())
    }
}

const HANDLED: &[&str] = &[
    "fastq_convert",
    "fasta2fastq",
    "sff_convert",
    "udb2fasta",
    "fastq_chars",
    "fastq_stats",
    "fastq_eestats",
    "fastq_eestats2",
    "fastx_filter",
    "fastq_filter",
    "fastq_mergepairs",
    "fastq_join",
    "fastx_syncpairs",
    "fastx_uniques",
    "derep_fulllength",
    "derep_id",
    "derep_prefix",
    "derep_smallmem",
    "rereplicate",
    "fastx_revcomp",
    "fastx_mask",
    "maskfasta",
    "cut",
    "scramble",
    "fastx_getseq",
    "fastx_getseqs",
    "fastx_getsubseq",
    "sortbylength",
    "sortbysize",
    "shuffle",
    "fastx_subsample",
    "usearch_global",
    "search_exact",
    "allpairs_global",
    "orient",
    "cluster_fast",
    "cluster_size",
    "cluster_smallmem",
    "cluster_unoise",
    "uchime_denovo",
    "uchime2_denovo",
    "uchime3_denovo",
    "uchime_ref",
    "chimeras_denovo",
    "sintax",
    "makeudb_usearch",
    "udbinfo",
    "udbstats",
];

/// The hit files of a search, made as they are first asked for: the
/// search with the hit choice given (`all`), with `--maxhits 1` (`top`),
/// or with neither (`every`).
struct Hits {
    tool: &'static str,
    queries: String,
    base: Vec<String>,
    choice: Vec<String>,
    cache: Vec<(Vec<String>, String)>,
    fixed: Option<String>,
}

impl Hits {
    fn new(tool: &'static str, queries: &str, base: Vec<String>, choice: Vec<String>) -> Hits {
        Hits {
            tool,
            queries: queries.to_string(),
            base,
            choice,
            cache: vec![],
            fixed: None,
        }
    }

    fn fixed(path: String) -> Hits {
        Hits {
            tool: "",
            queries: String::new(),
            base: vec![],
            choice: vec![],
            cache: vec![],
            fixed: Some(path),
        }
    }

    fn get(&mut self, t: &mut Translator, kind: &str) -> String {
        if let Some(f) = &self.fixed {
            return f.clone();
        }
        let key = match kind {
            "all" => self.choice.clone(),
            "top" => s(&["--maxhits", "1"]),
            _ => vec![],
        };
        if let Some((_, path)) = self.cache.iter().find(|(k, _)| *k == key) {
            return path.clone();
        }
        let mut args = self.base.clone();
        args.extend(key.iter().cloned());
        let path = t.step(self.tool, &self.queries.clone(), args, None, "uc");
        self.cache.push((key, path.clone()));
        path
    }
}

/// The rsearch pipeline for vsearch's command line `argv` (with or without
/// the program name), writing its files into `workdir`.
pub fn translate<S: AsRef<str>>(argv: &[S], workdir: &str) -> Result<Pipeline, TranslationError> {
    let a = Args::parse(argv)?;
    // vsearch's thread count, for rsearch's threaded tools; both default
    // to one thread per core
    let threads = a.peek("threads").and_then(|t| t.parse::<usize>().ok());
    let t = Translator {
        cmd: a.command.clone(),
        inp: a.input.clone().unwrap_or_default(),
        a,
        work: workdir.to_string(),
        p: Pipeline {
            command: String::new(),
            argv: argv.iter().map(|x| x.as_ref().to_string()).collect(),
            ..Default::default()
        },
        n: 0,
        width: Num {
            v: 80.0,
            text: "80".into(),
        },
        formats: HashMap::new(),
    };
    let mut t = t;
    t.p.command = t.cmd.clone();
    t.translate().map(|mut p| {
        if let Some(n) = threads {
            p.set_threads(n);
        }
        p
    })
}
