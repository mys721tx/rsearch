//! Runs the cases of sources over seeds on several threads, recording each
//! case as a line of JSON, and summarises recorded runs.
//!
//! The runner's threads only run cases. A reporter thread of its own takes
//! each seed's records from them, writes each record whole, and prints the
//! differences, so that what an interruption cuts short never reaches the
//! results.

use std::collections::{BTreeMap, VecDeque};
use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
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
    /// the threads of rsearch's threaded tools, whatever vsearch's argv says
    pub threads: Option<usize>,
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
        if STOP.load(Ordering::SeqCst) {
            break;
        }
        let t = Instant::now();
        let mut rec = run_case(
            case,
            &dir.join(format!("c{k}")),
            &opts.vsearch,
            &opts.bindir,
            opts.keep,
            opts.threads,
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
    let file = match &opts.out {
        Some(p) => {
            let (file, held) = open_out(p)?;
            if held > 0 {
                eprintln!(
                    "Warning: {} already holds {held} records; this run's are added after them",
                    p.display()
                );
            }
            Some(file)
        }
        None => None,
    };
    let queue = Arc::new(Mutex::new(items));
    let opts = Arc::new(opts);
    let (tx, rx) = mpsc::channel();
    let reporter = std::thread::spawn(move || report(rx, file, &STOP));
    let mut handles = Vec::new();
    for _ in 0..opts.jobs.max(1) {
        let (queue, opts, tx) = (queue.clone(), opts.clone(), tx.clone());
        handles.push(std::thread::spawn(move || {
            while !STOP.load(Ordering::SeqCst) {
                let item = queue.lock().unwrap().pop_front();
                let Some((name, seed)) = item else { break };
                if tx.send(work(&name, seed, &opts)).is_err() {
                    break;
                }
            }
        }));
    }
    drop(tx);
    for h in handles {
        let _ = h.join();
    }
    reporter.join().expect("the reporter thread")
}

/// Set when the run is interrupted: no case starts after it, and no record
/// that arrives after it is written.
static STOP: AtomicBool = AtomicBool::new(false);

/// whether the run was interrupted.
pub fn interrupted() -> bool {
    STOP.load(Ordering::SeqCst)
}

/// makes SIGINT and SIGTERM stop the run instead of the process, so that
/// the reporter ends the results on a whole record.
pub fn catch_interrupts() {
    #[cfg(unix)]
    {
        extern "C" fn on_signal(_: libc::c_int) {
            // only an atomic store: all a signal handler may safely do
            STOP.store(true, Ordering::SeqCst);
        }
        let handler = on_signal as extern "C" fn(libc::c_int) as libc::sighandler_t;
        // SAFETY: the handler only stores to an atomic
        unsafe {
            libc::signal(libc::SIGINT, handler);
            libc::signal(libc::SIGTERM, handler);
        }
    }
}

/// The reporter: takes each seed's records from `rx` until the runner's
/// threads are done, counts them, writes each whole to `file`, and prints
/// the first 200 differences. Records that arrive once `stop` is set come
/// from cases the interruption cut short, and are dropped.
fn report(
    rx: mpsc::Receiver<Vec<Map<String, Value>>>,
    mut file: Option<std::fs::File>,
    stop: &AtomicBool,
) -> std::io::Result<Tally> {
    let mut tally = Tally::new();
    let mut shown = 0;
    for recs in rx {
        if stop.load(Ordering::SeqCst) {
            continue;
        }
        for r in recs {
            let source = r["source"].as_str().unwrap_or("").to_string();
            let outcome = r["outcome"].as_str().unwrap_or("").to_string();
            *tally
                .entry(source.clone())
                .or_default()
                .entry(outcome.clone())
                .or_default() += 1;
            if let Some(f) = file.as_mut() {
                write_record(f, &r)?;
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
    Ok(tally)
}

/// opens `path` to add records to, creating it, and returns it with the
/// number of records it already holds. When it ends partway through a line,
/// as a run stopped while writing leaves it, a newline goes first, so that
/// the next record starts a line of its own.
fn open_out(path: &Path) -> std::io::Result<(std::fs::File, usize)> {
    let mut file = OpenOptions::new()
        .create(true)
        .read(true)
        .append(true)
        .open(path)?;
    let mut held = 0;
    let mut last = b'\n';
    let mut reader = BufReader::new(&file);
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            break;
        }
        held += buffer.iter().filter(|&&b| b == b'\n').count();
        last = *buffer.last().unwrap();
        let n = buffer.len();
        reader.consume(n);
    }
    if last != b'\n' {
        held += 1;
        file.write_all(b"\n")?;
    }
    Ok((file, held))
}

/// writes `record` as a line of JSON, in one write.
fn write_record(file: &mut std::fs::File, record: &Map<String, Value>) -> std::io::Result<()> {
    let mut line = Value::Object(record.clone()).to_string();
    line.push('\n');
    file.write_all(line.as_bytes())
}

/// the records of `paths`, each case (source, seed, index) once, by its
/// last record, and warnings for the lines that are not records and for
/// the cases recorded more than once.
fn read_records(paths: &[PathBuf]) -> std::io::Result<(Vec<Value>, Vec<String>)> {
    let mut records: Vec<Value> = Vec::new();
    let mut place: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut warnings = Vec::new();
    let mut repeated = 0;
    for path in paths {
        for (n, line) in BufReader::new(std::fs::File::open(path)?)
            .lines()
            .enumerate()
        {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let record: Value = match serde_json::from_str(&line) {
                Ok(r) => r,
                Err(e) => {
                    warnings.push(format!(
                        "{}:{}: not a record, skipped ({e})",
                        path.display(),
                        n + 1
                    ));
                    continue;
                }
            };
            let key = format!(
                "{} {} {}",
                record["source"], record["seed"], record["index"]
            );
            match place.get(&key) {
                Some(&at) => {
                    if records[at] != Value::Null {
                        repeated += 1;
                    }
                    records[at] = Value::Null;
                    place.insert(key, records.len());
                    records.push(record);
                }
                None => {
                    place.insert(key, records.len());
                    records.push(record);
                }
            }
        }
    }
    if repeated > 0 {
        let cases = if repeated == 1 {
            "case is"
        } else {
            "cases are"
        };
        warnings.push(format!(
            "{repeated} {cases} recorded more than once; the last record of each is counted"
        ));
    }
    records.retain(|r| *r != Value::Null);
    Ok((records, warnings))
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
    let (records, warnings) = read_records(paths)?;
    for warning in &warnings {
        eprintln!("Warning: {warning}");
    }
    for r in records {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn temp(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("rscompare-runner-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(name)
    }

    fn record(source: &str, seed: u64, index: u64, outcome: &str) -> Map<String, Value> {
        let Value::Object(m) = json!({
            "case": "c", "source": source, "seed": seed, "index": index, "outcome": outcome,
        }) else {
            unreachable!()
        };
        m
    }

    #[test]
    fn reports_whole_batches_until_stopped() {
        // the reporter writes each batch whole, then drops what arrives
        // after a stop: cases cut short by the interruption are not results
        let path = temp("reporter.jsonl");
        let _ = std::fs::remove_file(&path);
        let (file, _) = open_out(&path).unwrap();
        let stop = std::sync::atomic::AtomicBool::new(false);
        let (tx, rx) = mpsc::channel();
        tx.send(vec![
            record("derep", 1, 0, "pass"),
            record("derep", 1, 1, "fail"),
        ])
        .unwrap();
        tx.send(vec![record("conv", 2, 0, "pass")]).unwrap();
        let sender = std::thread::spawn(move || tx);
        let tx = sender.join().unwrap();
        std::thread::scope(|scope| {
            let reporter = scope.spawn(|| report(rx, Some(file), &stop));
            // let the reporter take the first two batches
            while std::fs::read_to_string(&path).unwrap().lines().count() < 3 {
                std::thread::yield_now();
            }
            stop.store(true, std::sync::atomic::Ordering::SeqCst);
            tx.send(vec![record("relab", 3, 0, "fail")]).unwrap();
            drop(tx);
            let tally = reporter.join().unwrap().unwrap();
            assert_eq!(tally["derep"]["pass"], 1);
            assert_eq!(tally["derep"]["fail"], 1);
            assert_eq!(tally["conv"]["pass"], 1);
            assert!(!tally.contains_key("relab"));
        });
        let text = std::fs::read_to_string(&path).unwrap();
        assert_eq!(text.lines().count(), 3);
        assert!(text.ends_with('\n'));
        assert!(
            text.lines()
                .all(|l| serde_json::from_str::<Value>(l).is_ok())
        );
    }

    #[test]
    fn starts_a_line_after_a_record_cut_short() {
        // a run stopped while writing leaves half a line; the next run's
        // first record must not join it
        let path = temp("cut.jsonl");
        std::fs::write(&path, "{\"case\":\"a\"}\n{\"case\":\"b\",\"sou").unwrap();
        let (mut file, held) = open_out(&path).unwrap();
        assert_eq!(held, 2);
        write_record(&mut file, &record("derep", 1, 0, "pass")).unwrap();
        drop(file);
        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[1], "{\"case\":\"b\",\"sou");
        assert!(lines[2].starts_with('{') && serde_json::from_str::<Value>(lines[2]).is_ok());
        assert!(text.ends_with('\n'));
    }

    #[test]
    fn creates_and_counts_what_out_holds() {
        let path = temp("new.jsonl");
        let _ = std::fs::remove_file(&path);
        let (mut file, held) = open_out(&path).unwrap();
        assert_eq!(held, 0);
        for k in 0..3 {
            write_record(&mut file, &record("conv", 2, k, "pass")).unwrap();
        }
        drop(file);
        assert_eq!(open_out(&path).unwrap().1, 3);
        assert_eq!(std::fs::read_to_string(&path).unwrap().lines().count(), 3);
    }

    #[test]
    fn reads_records_skipping_lines_that_are_not_and_counting_each_case_once() {
        let path = temp("mixed.jsonl");
        let line = |m: Map<String, Value>| format!("{}\n", Value::Object(m));
        let text = [
            line(record("derep", 1, 0, "fail")),
            line(record("derep", 1, 1, "pass")),
            "{\"case\":\"cut\",\"source\":\"uch{\"case\":\"x\"}\n".to_string(),
            "\n".to_string(),
            line(record("derep", 1, 0, "pass")),
            line(record("conv", 3, 0, "skip")),
        ]
        .concat();
        std::fs::write(&path, text).unwrap();
        let (records, warnings) = read_records(std::slice::from_ref(&path)).unwrap();
        // the last record of derep 1 #0 counts, once
        let outcomes: Vec<(String, String)> = records
            .iter()
            .map(|r| {
                (
                    format!("{} {} {}", r["source"], r["seed"], r["index"]),
                    r["outcome"].as_str().unwrap().to_string(),
                )
            })
            .collect();
        assert_eq!(outcomes.len(), 3);
        assert!(outcomes.contains(&("\"derep\" 1 0".into(), "pass".into())));
        assert!(outcomes.contains(&("\"derep\" 1 1".into(), "pass".into())));
        assert!(outcomes.contains(&("\"conv\" 3 0".into(), "skip".into())));
        assert_eq!(warnings.len(), 2, "{warnings:?}");
        assert!(
            warnings[0].starts_with(&format!("{}:3: not a record, skipped", path.display())),
            "{warnings:?}"
        );
        assert_eq!(
            warnings[1],
            "1 case is recorded more than once; the last record of each is counted"
        );
    }
}
