//! Runs the cases of sources over seeds on several threads, recording each
//! case as a line of JSON, and summarises recorded runs.

use std::collections::{BTreeMap, VecDeque};
use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use serde_json::{Map, Value, json};

use crate::cases::run_case;
use crate::sources::{self, Env};

/// Outcome counts per source.
pub type Tally = BTreeMap<String, BTreeMap<String, u64>>;

/// What a run needs.
pub struct Options {
    pub root: PathBuf,
    pub vsearch: PathBuf,
    pub bindir: PathBuf,
    pub out: Option<PathBuf>,
    pub keep: bool,
    pub jobs: usize,
}

fn work(name: &str, seed: u64, opts: &Options) -> Vec<Map<String, Value>> {
    let src = sources::find(name).expect("a known source");
    let dir = opts.root.join(format!("{name}-{seed}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let base = |mut m: Map<String, Value>| {
        m.insert("source".into(), json!(name));
        m.insert("seed".into(), json!(seed));
        m
    };
    let env = Env {
        bindir: opts.bindir.clone(),
        vsearch: opts.vsearch.clone(),
    };
    let cases = match (src.make)(&dir, seed, &env) {
        Ok(c) => c,
        Err(e) => {
            let mut m = Map::new();
            m.insert("case".into(), json!("(generator)"));
            m.insert("outcome".into(), json!("error"));
            m.insert("reason".into(), json!(e.to_string()));
            return vec![base(m)];
        }
    };
    let mut out = Vec::new();
    for (k, case) in cases.iter().enumerate() {
        let t = Instant::now();
        let mut rec = run_case(
            case,
            &dir.join(format!("c{k}")),
            &opts.vsearch,
            &opts.bindir,
            opts.keep,
        );
        rec.insert("index".into(), json!(k));
        rec.insert(
            "wall".into(),
            json!((t.elapsed().as_secs_f64() * 1e3).round() / 1e3),
        );
        out.push(base(rec));
    }
    if !opts.keep
        && out
            .iter()
            .all(|r| matches!(r["outcome"].as_str(), Some("pass" | "skip")))
    {
        let _ = std::fs::remove_dir_all(&dir);
    }
    out
}

/// Runs the sources' cases at the seeds (each source's own by default),
/// printing a DIFF line for each failure (the first 200), and returns the
/// tally.
pub fn run(names: &[String], seeds: Option<&[u64]>, opts: Options) -> std::io::Result<Tally> {
    std::fs::create_dir_all(&opts.root)?;
    let mut items = VecDeque::new();
    for name in names {
        let src = sources::find(name).expect("a known source");
        let list = seeds.map_or_else(|| (src.seeds)(), |s| s.to_vec());
        for seed in list {
            items.push_back((name.clone(), seed));
        }
    }
    let queue = Arc::new(Mutex::new(items));
    let opts = Arc::new(opts);
    let (tx, rx) = mpsc::channel();
    let mut handles = Vec::new();
    for _ in 0..opts.jobs.max(1) {
        let (queue, opts, tx) = (queue.clone(), opts.clone(), tx.clone());
        handles.push(std::thread::spawn(move || {
            loop {
                let item = queue.lock().unwrap().pop_front();
                let Some((name, seed)) = item else { break };
                if tx.send(work(&name, seed, &opts)).is_err() {
                    break;
                }
            }
        }));
    }
    drop(tx);
    let mut file = match &opts.out {
        Some(p) => Some(OpenOptions::new().create(true).append(true).open(p)?),
        None => None,
    };
    let mut tally = Tally::new();
    let mut shown = 0;
    for recs in rx {
        for r in recs {
            let source = r["source"].as_str().unwrap_or("").to_string();
            let outcome = r["outcome"].as_str().unwrap_or("").to_string();
            *tally
                .entry(source.clone())
                .or_default()
                .entry(outcome.clone())
                .or_default() += 1;
            if let Some(f) = file.as_mut() {
                writeln!(f, "{}", Value::Object(r.clone()))?;
            }
            if (outcome == "fail" || outcome == "error") && shown < 200 {
                shown += 1;
                let what: Vec<String> = match r.get("diffs") {
                    Some(Value::Array(d)) => d
                        .iter()
                        .map(|x| x.as_str().unwrap_or("").to_string())
                        .collect(),
                    _ => vec![
                        r.get("reason")
                            .and_then(Value::as_str)
                            .unwrap_or("")
                            .to_string(),
                    ],
                };
                let text: String = what.join("; ").chars().take(400).collect();
                println!(
                    "DIFF {source} seed {} [{}]: {text}",
                    r["seed"],
                    r["case"].as_str().unwrap_or("")
                );
            }
        }
        if let Some(f) = file.as_mut() {
            f.flush()?;
        }
    }
    for h in handles {
        let _ = h.join();
    }
    Ok(tally)
}

const OUTCOMES: [&str; 4] = ["pass", "fail", "skip", "error"];

/// Prints the outcome counts per source, and their total.
pub fn print_tally(tally: &Tally) {
    let mut total: BTreeMap<String, u64> = BTreeMap::new();
    let line = |name: &str, t: &BTreeMap<String, u64>| {
        let parts: Vec<String> = OUTCOMES
            .iter()
            .filter_map(|k| t.get(*k).map(|n| format!("{k}={n}")))
            .collect();
        println!("{name:12} {}", parts.join(" "));
    };
    for (name, t) in tally {
        for (k, v) in t {
            *total.entry(k.clone()).or_default() += v;
        }
        line(name, t);
    }
    if tally.len() > 1 {
        line("total", &total);
    }
}

/// Outcomes, coverage tallies, times and caveats of recorded runs.
pub fn summary(paths: &[PathBuf]) -> std::io::Result<()> {
    let mut tally = Tally::new();
    let mut cover: BTreeMap<String, BTreeMap<String, BTreeMap<String, u64>>> = BTreeMap::new();
    let mut caveats: BTreeMap<String, u64> = BTreeMap::new();
    let mut seconds: BTreeMap<String, (f64, f64)> = BTreeMap::new();
    for path in paths {
        for line in BufReader::new(std::fs::File::open(path)?).lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let r: Value = serde_json::from_str(&line).map_err(std::io::Error::other)?;
            let s = r["source"].as_str().unwrap_or("").to_string();
            let outcome = r["outcome"].as_str().unwrap_or("").to_string();
            *tally
                .entry(s.clone())
                .or_default()
                .entry(outcome)
                .or_default() += 1;
            let c = cover.entry(s.clone()).or_default();
            if let Some(Value::Object(cov)) = r.get("coverage") {
                for (k, v) in cov {
                    let v = match v {
                        Value::Object(o) => o
                            .iter()
                            .map(|(a, b)| format!("{a}{b}"))
                            .collect::<Vec<_>>()
                            .join(","),
                        Value::Number(n) if k != "exit" => match n.as_f64() {
                            Some(0.0) => "0".into(),
                            Some(1.0) => "1".into(),
                            _ => ">1".into(),
                        },
                        Value::String(t) => t.clone(),
                        other => other.to_string(),
                    };
                    *c.entry(k.clone()).or_default().entry(v).or_default() += 1;
                }
            }
            if let Some(Value::Array(cs)) = r.get("caveats") {
                for cv in cs {
                    *caveats
                        .entry(cv.as_str().unwrap_or("").to_string())
                        .or_default() += 1;
                }
            }
            if let Some(sec) = r.get("seconds") {
                let e = seconds.entry(s).or_default();
                e.0 += sec["vsearch"].as_f64().unwrap_or(0.0);
                e.1 += sec["rsearch"].as_f64().unwrap_or(0.0);
            }
        }
    }
    print_tally(&tally);
    println!();
    for (s, keys) in &cover {
        println!("{s}:");
        for (k, vals) in keys {
            let mut top: Vec<(&String, &u64)> = vals.iter().collect();
            top.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
            let parts: Vec<String> = top
                .iter()
                .take(8)
                .map(|(v, n)| format!("{v}={n}"))
                .collect();
            println!("  {k}: {}", parts.join(", "));
        }
        if let Some((v, r)) = seconds.get(s) {
            println!("  time: vsearch {v:.1}s, rsearch {r:.1}s");
        }
    }
    if !caveats.is_empty() {
        println!("\ncaveats:");
        let mut list: Vec<(&String, &u64)> = caveats.iter().collect();
        list.sort_by(|a, b| b.1.cmp(a.1));
        for (cv, n) in list {
            println!("  {n:7}  {cv}");
        }
    }
    Ok(())
}

/// Seeds from text such as `1-50,70`.
pub fn parse_seeds(text: &str) -> Result<Vec<u64>, String> {
    let mut out = Vec::new();
    for part in text.split(',').filter(|p| !p.is_empty()) {
        match part.split_once('-') {
            Some((a, b)) => {
                let a: u64 = a.parse().map_err(|_| format!("bad seed range {part}"))?;
                let b: u64 = b.parse().map_err(|_| format!("bad seed range {part}"))?;
                out.extend(a..=b);
            }
            None => out.push(part.parse().map_err(|_| format!("bad seed {part}"))?),
        }
    }
    Ok(out)
}

/// The number of threads to use by default.
pub fn default_jobs() -> usize {
    std::thread::available_parallelism().map_or(4, |n| n.get())
}

/// The default directory for case files.
pub fn default_root() -> PathBuf {
    std::env::var_os("TMPDIR")
        .map_or_else(|| PathBuf::from("/tmp"), PathBuf::from)
        .join("rscompare")
}

/// The repository root (the directory above `compare/`).
pub fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}
