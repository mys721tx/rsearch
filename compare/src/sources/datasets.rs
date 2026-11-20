//! Real datasets: the manifest `compare/datasets.json`, downloaded by
//! `rscompare fetch` into a data directory (`$RSCOMPARE_DATA`,
//! by default `$TMPDIR/rsearch-data`, one directory per dataset).
//!
//! The source reads each manifest entry's vsearch command lines and runs
//! them as cases; the seed is the entry's place in the manifest (from 1).
//! An entry is a JSON object; the source uses:
//!
//! - its name: `id`, else `name`, else `accession`;
//! - its command lines: `vsearch` or `commands` (or `vsearch_commands`), a
//!   list of strings or of argv lists, without the program name. In them
//!   `{dir}` is the dataset's directory (relative file names are taken
//!   from it too); `{O}` is the output directory. Output options without
//!   `{O}` are moved there, named by the option.
//!
//! Every command line gets `--threads 1` and `--quiet`, and
//! `--notrunclabels` where vsearch would cut labels. An entry whose files
//! are missing is skipped, naming them.

use std::io;
use std::path::{Path, PathBuf};

use serde_json::Value;

use super::{Env, Source, shell_split};
use crate::cases::{Case, OUT};
use crate::options::Args;
use crate::translate::WHOLE_HEADERS;

pub fn sources() -> Vec<Source> {
    vec![Source {
        name: "datasets",
        seeds: || (1..=entries().len() as u64).collect(),
        about: "the real datasets of compare/datasets.json, by their place in it (seeds 1 to N)",
        make: datasets,
    }]
}

const OUTPUTS: &[&str] = &[
    "alnout",
    "biomout",
    "blast6out",
    "borderline",
    "centroids",
    "chimeras",
    "clusters",
    "consout",
    "dbmatched",
    "dbnotmatched",
    "eetabbedout",
    "fastaout",
    "fastaout_discarded",
    "fastaout_discarded_rev",
    "fastaout_notmerged_fwd",
    "fastaout_notmerged_rev",
    "fastaout_orphans",
    "fastaout_orphans_rev",
    "fastaout_rev",
    "fastapairs",
    "fastqout",
    "fastqout_discarded",
    "fastqout_discarded_rev",
    "fastqout_notmerged_fwd",
    "fastqout_notmerged_rev",
    "fastqout_orphans",
    "fastqout_orphans_rev",
    "fastqout_rev",
    "lcaout",
    "log",
    "matched",
    "mothur_shared_out",
    "msaout",
    "nonchimeras",
    "notmatched",
    "notmatchedfq",
    "otutabout",
    "output",
    "profile",
    "qsegout",
    "samout",
    "tabbedout",
    "tsegout",
    "uc",
    "uchimealns",
    "uchimeout",
    "userout",
];

/// Options whose value names an input file.
const INPUTS: &[&str] = &["db", "reverse", "labels", "label_words"];

/// The manifest.
pub fn manifest() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("datasets.json")
}

/// The manifest's entries, or none.
pub fn entries() -> Vec<Value> {
    let Ok(text) = std::fs::read_to_string(manifest()) else {
        return vec![];
    };
    let Ok(data) = serde_json::from_str::<Value>(&text) else {
        return vec![];
    };
    let list = match data {
        Value::Array(a) => a,
        Value::Object(mut o) => match o.remove("datasets") {
            Some(Value::Array(a)) => a,
            _ => o.into_iter().map(|(_, v)| v).collect(),
        },
        _ => vec![],
    };
    list.into_iter().filter(Value::is_object).collect()
}

fn data_dir() -> PathBuf {
    std::env::var_os("RSCOMPARE_DATA")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            std::env::var_os("TMPDIR")
                .map_or_else(|| PathBuf::from("/tmp"), PathBuf::from)
                .join("rsearch-data")
        })
}

fn entry_name(e: &Value) -> String {
    ["id", "name", "accession"]
        .iter()
        .find_map(|k| e.get(*k).filter(|v| !v.is_null()))
        .map_or("?".into(), |v| {
            v.as_str().map_or_else(|| v.to_string(), str::to_string)
        })
}

/// A manifest command line as a case's argv: inputs from the dataset's
/// directory, outputs moved under [`OUT`]; and the inputs missing.
pub fn prepare(cmd: &Value, dir: &str) -> Result<(Vec<String>, Vec<String>), String> {
    let mut words: Vec<String> = match cmd {
        Value::String(s) => shell_split(&s.replace("{dir}", dir)),
        Value::Array(a) => a
            .iter()
            .map(|w| {
                w.as_str()
                    .map_or_else(|| w.to_string(), str::to_string)
                    .replace("{dir}", dir)
            })
            .collect(),
        _ => return Err("not a command line".into()),
    };
    if words.first().is_some_and(|w| !w.starts_with('-')) {
        words.remove(0); // the program
    }
    let args = Args::parse(&words).map_err(|e| e.to_string())?;
    let (mut out, mut missing) = (Vec::new(), Vec::new());
    let mut k = 0;
    while k < words.len() {
        let w = &words[k];
        let name = w.trim_start_matches('-');
        out.push(w.clone());
        if k + 1 < words.len() {
            let v = &words[k + 1];
            if OUTPUTS.contains(&name) {
                out.push(if v.contains("{O}") {
                    v.replace("{O}", OUT)
                } else {
                    format!("{OUT}/{name}")
                });
                k += 2;
                continue;
            }
            if name == args.command || INPUTS.contains(&name) {
                let path = if Path::new(v).is_absolute() {
                    v.clone()
                } else {
                    format!("{dir}/{v}")
                };
                if !Path::new(&path).exists() {
                    missing.push(path.clone());
                }
                out.push(path);
                k += 2;
                continue;
            }
        }
        k += 1;
    }
    let has = |out: &[String], n: &str| out.iter().any(|w| w.trim_start_matches('-') == n);
    if !has(&out, "threads") {
        out.extend(["--threads".to_string(), "1".to_string()]);
    }
    if !has(&out, "quiet") {
        out.push("--quiet".into());
    }
    if !WHOLE_HEADERS.contains(&args.command.as_str()) && !has(&out, "notrunclabels") {
        out.push("--notrunclabels".into());
    }
    Ok((out, missing))
}

fn commands(e: &Value) -> Vec<Value> {
    for key in ["vsearch", "commands", "vsearch_commands"] {
        match e.get(key) {
            Some(Value::Array(a)) if !a.is_empty() => return a.clone(),
            Some(v @ Value::String(_)) => return vec![v.clone()],
            _ => {}
        }
    }
    vec![]
}

fn datasets(_: &Path, seed: u64, _: &Env) -> io::Result<Vec<Case>> {
    let es = entries();
    if seed == 0 || seed as usize > es.len() {
        return Ok(vec![]);
    }
    let e = &es[seed as usize - 1];
    let name = entry_name(e);
    let dir = data_dir().join(&name).to_string_lossy().into_owned();
    let mut cases = Vec::new();
    for (k, cmd) in commands(e).iter().enumerate() {
        let label = format!("{name} #{k}");
        match prepare(cmd, &dir) {
            Ok((argv, missing)) if missing.is_empty() => cases.push(
                Case::new(label, argv)
                    .errors("exact")
                    .warn(&["stripped"])
                    .cover("dataset", name.as_str()),
            ),
            Ok((_, missing)) => {
                let why = format!(
                    "missing {}",
                    missing
                        .iter()
                        .take(3)
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(", ")
                );
                cases.push(Case::new(label, vec![]).skip(Box::new(move |_, _| Some(why.clone()))));
            }
            Err(e) => {
                let why = format!("unparsable command line: {e}");
                cases.push(Case::new(label, vec![]).skip(Box::new(move |_, _| Some(why.clone()))));
            }
        }
    }
    Ok(cases)
}
