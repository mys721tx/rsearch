//! A case, and running one. vsearch runs on its command line, and rsearch
//! on the translated pipeline (or, for a probe, on the steps the case
//! gives). Their outputs, exit codes, errors and warnings are then
//! compared.

use std::collections::HashMap;
use std::path::Path;

use serde_json::{Map, Value, json};

use crate::comparators as C;
use crate::execute::{Run, run_steps, run_vsearch};
use crate::translate::{Output, Pipeline, Step, translate};

/// Stands for the case's vsearch output directory in its argv.
pub const OUT: &str = "@OUT@";

/// A probe: rsearch steps the case gives instead of the translation, made
/// for the rsearch directory, with the file standing for each vsearch
/// output option.
pub type Probe = Box<dyn Fn(&Path) -> (Vec<Step>, Vec<(String, String)>)>;
/// A rule that skips a case after both ran, with the reason.
pub type Skip = Box<dyn Fn(Option<&Run>, Option<&Run>) -> Option<String>>;

/// One comparison.
pub struct Case {
    pub name: String,
    /// vsearch's argv, outputs under [`OUT`].
    pub vsearch: Vec<String>,
    /// `exit`, or an error kind of [`comparators::error_kinds`](crate::comparators::error_kinds).
    pub errors: String,
    pub warnings: Vec<String>,
    /// Compare outputs when vsearch fails.
    pub outputs_on_error: bool,
    /// (vsearch's exit, rsearch's exit) taken as equal to the value.
    pub exit_map: Vec<((i32, i32), i32)>,
    pub probe: Option<Probe>,
    /// The comparator per vsearch option, overriding the translation's.
    pub comparators: HashMap<String, String>,
    /// Outputs and warnings compared only when vsearch succeeds.
    pub success_only: Vec<String>,
    pub skip: Option<Skip>,
    pub coverage: Map<String, Value>,
    /// Counts from vsearch's outputs for the coverage record: name, the
    /// option, and `records`, `uc`, `lines` or `col N VALUE`.
    pub counts: Vec<(String, String, String)>,
}

impl Case {
    pub fn new(name: impl Into<String>, vsearch: Vec<String>) -> Case {
        Case {
            name: name.into(),
            vsearch,
            errors: "exit".into(),
            warnings: vec![],
            outputs_on_error: true,
            exit_map: vec![],
            probe: None,
            comparators: HashMap::new(),
            success_only: vec![],
            skip: None,
            coverage: Map::new(),
            counts: vec![],
        }
    }

    pub fn errors(mut self, kind: &str) -> Case {
        self.errors = kind.into();
        self
    }

    pub fn warn(mut self, names: &[&str]) -> Case {
        self.warnings = names.iter().map(|s| s.to_string()).collect();
        self
    }

    pub fn success_outputs_only(mut self) -> Case {
        self.outputs_on_error = false;
        self
    }

    pub fn exit_map(mut self, v: i32, r: i32, as_: i32) -> Case {
        self.exit_map.push(((v, r), as_));
        self
    }

    pub fn probe(mut self, p: Probe) -> Case {
        self.probe = Some(p);
        self
    }

    pub fn comparator(mut self, option: &str, name: &str) -> Case {
        self.comparators.insert(option.into(), name.into());
        self
    }

    pub fn success_only(mut self, names: &[&str]) -> Case {
        self.success_only = names.iter().map(|s| s.to_string()).collect();
        self
    }

    pub fn skip(mut self, s: Skip) -> Case {
        self.skip = Some(s);
        self
    }

    pub fn cover(mut self, key: &str, value: impl Into<Value>) -> Case {
        self.coverage.insert(key.into(), value.into());
        self
    }

    pub fn count(mut self, name: &str, option: &str, kind: &str) -> Case {
        self.counts.push((name.into(), option.into(), kind.into()));
        self
    }
}

fn subst(argv: &[String], vdir: &str) -> Vec<String> {
    argv.iter().map(|a| a.replace(OUT, vdir)).collect()
}

fn option_value<'a>(argv: &'a [String], option: &str) -> Option<&'a str> {
    argv.windows(2)
        .find(|w| w[0].trim_start_matches('-') == option)
        .map(|w| w[1].as_str())
}

fn count(path: &str, kind: &str) -> Value {
    let Ok(data) = std::fs::read(path) else {
        return Value::Null;
    };
    let mut lines: Vec<&[u8]> = data.split(|&b| b == b'\n').collect();
    if lines.last() == Some(&&b""[..]) {
        lines.pop();
    }
    match kind {
        "lines" => json!(lines.len()),
        "records" if data.is_empty() => json!(0),
        "records" if data[0] == b'@' => json!(lines.len() / 4),
        "records" => json!(lines.iter().filter(|l| l.starts_with(b">")).count()),
        "uc" => {
            let mut out = Map::new();
            for l in lines {
                let key = String::from_utf8_lossy(&l[..l.len().min(1)]).into_owned();
                let n = out.get(&key).and_then(Value::as_u64).unwrap_or(0);
                out.insert(key, json!(n + 1));
            }
            Value::Object(out)
        }
        k if k.starts_with("col ") => {
            let mut it = k.splitn(3, ' ').skip(1);
            let n: usize = it.next().and_then(|n| n.parse().ok()).unwrap_or(0);
            let value = it.next().unwrap_or("").as_bytes();
            json!(
                lines
                    .iter()
                    .filter(|l| l.split(|&b| b == b'\t').nth(n) == Some(value))
                    .count()
            )
        }
        _ => Value::Null,
    }
}

fn round4(x: f64) -> f64 {
    (x * 1e4).round() / 1e4
}

/// Runs a case in `workdir` (its `v/` holds vsearch's outputs, `r/`
/// rsearch's) and returns its record. The directory is removed when the
/// case passes, unless `keep`.
pub fn run_case(
    case: &Case,
    workdir: &Path,
    vsearch: &Path,
    bindir: &Path,
    keep: bool,
    threads: Option<usize>,
) -> Map<String, Value> {
    let mut rec = Map::new();
    rec.insert("case".into(), json!(case.name));
    if case.vsearch.is_empty() {
        // nothing to run (a dataset not downloaded)
        let why = case
            .skip
            .as_ref()
            .and_then(|s| s(None, None))
            .unwrap_or("no command".into());
        rec.insert("outcome".into(), json!("skip"));
        rec.insert("reason".into(), json!(why));
        return rec;
    }
    let (vdir, rdir) = (workdir.join("v"), workdir.join("r"));
    let _ = std::fs::remove_dir_all(workdir);
    std::fs::create_dir_all(&vdir).unwrap();
    std::fs::create_dir_all(&rdir).unwrap();
    let argv = subst(&case.vsearch, &vdir.to_string_lossy());
    rec.insert("vsearch".into(), json!(argv));
    let mut coverage = case.coverage.clone();
    let pipe = if let Some(probe) = &case.probe {
        let (steps, outmap) = probe(&rdir);
        rec.insert("kind".into(), json!("probe"));
        Pipeline {
            command: "probe".into(),
            argv: argv.clone(),
            steps,
            outputs: outmap
                .into_iter()
                .map(|(option, rsearch)| Output {
                    vsearch: option_value(&argv, &option).unwrap_or("").to_string(),
                    option,
                    rsearch,
                    comparator: "bytes".into(),
                })
                .collect(),
            caveats: vec![],
        }
    } else {
        match translate(&argv, &rdir.to_string_lossy()) {
            Ok(p) => {
                rec.insert("kind".into(), json!("translated"));
                p
            }
            Err(e) => {
                rec.insert("coverage".into(), Value::Object(coverage));
                rec.insert("outcome".into(), json!("error"));
                rec.insert("reason".into(), json!(format!("translation: {e}")));
                return rec;
            }
        }
    };
    let mut pipe = pipe;
    if let Some(n) = threads {
        pipe.set_threads(n);
    }
    rec.insert(
        "rsearch".into(),
        json!(pipe.render().split('\n').collect::<Vec<_>>()),
    );
    rec.insert("caveats".into(), json!(pipe.caveats));
    let vr = run_vsearch(vsearch, &argv, workdir);
    let rr = run_steps(&pipe.steps, bindir, workdir, &rdir);
    rec.insert(
        "exit".into(),
        json!({"vsearch": vr.exit, "rsearch": rr.exit}),
    );
    rec.insert(
        "seconds".into(),
        json!({"vsearch": round4(vr.seconds), "rsearch": round4(rr.seconds)}),
    );
    if let Some(f) = &rr.failed {
        rec.insert("failed_step".into(), json!(f));
    }
    for (name, opt, kind) in &case.counts {
        if let Some(path) = option_value(&argv, opt) {
            coverage.insert(name.clone(), count(path, kind));
        }
    }
    let verr = C::vsearch_error(&vr.stderr);
    coverage.insert("exit".into(), json!(vr.exit));
    if !verr.is_empty() {
        coverage.insert(
            "error".into(),
            json!(verr.chars().take(60).collect::<String>()),
        );
    }
    rec.insert("coverage".into(), Value::Object(coverage));
    if let Some(skip) = &case.skip
        && let Some(why) = skip(Some(&vr), Some(&rr))
    {
        rec.insert("outcome".into(), json!("skip"));
        rec.insert("reason".into(), json!(why));
        if !keep {
            let _ = std::fs::remove_dir_all(workdir);
        }
        return rec;
    }
    let mut diffs = Vec::new();
    let rexit = case
        .exit_map
        .iter()
        .find(|(k, _)| *k == (vr.exit, rr.exit))
        .map_or(rr.exit, |(_, v)| *v);
    if vr.exit != rexit {
        diffs.push(format!("exit: vsearch {}, rsearch {}", vr.exit, rr.exit));
    }
    if case.errors != "exit" {
        let (v, r) = C::error_kinds(&case.errors, &vr.stderr, &rr.stderr);
        if v != r {
            diffs.push(format!("error: vsearch {v:?} | rsearch {r:?}"));
        }
    }
    for w in &case.warnings {
        if vr.exit != 0 && case.success_only.contains(w) {
            continue;
        }
        let (v, r) = C::warnings(w, &vr.all_stderr, &rr.all_stderr);
        if v != r {
            diffs.push(format!("warning {w}: vsearch {v:?} | rsearch {r:?}"));
        }
    }
    if vr.exit == 0 || case.outputs_on_error {
        for o in &pipe.outputs {
            if vr.exit != 0 && case.success_only.contains(&o.option) {
                continue;
            }
            let comp = case.comparators.get(&o.option).unwrap_or(&o.comparator);
            if let Some(d) = C::compare_files(comp, &o.vsearch, &o.rsearch) {
                diffs.push(format!("--{} ({comp}): {d}", o.option));
            }
        }
    }
    rec.insert(
        "outcome".into(),
        json!(if diffs.is_empty() { "pass" } else { "fail" }),
    );
    if !diffs.is_empty() {
        rec.insert("diffs".into(), json!(diffs));
        rec.insert(
            "stderr".into(),
            json!({"vsearch": C::vsearch_error(&vr.stderr), "rsearch": C::rsearch_error(&rr.stderr)}),
        );
    }
    if !keep && diffs.is_empty() {
        let _ = std::fs::remove_dir_all(workdir);
    }
    rec
}
