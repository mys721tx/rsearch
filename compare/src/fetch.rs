//! Downloads the real datasets of `datasets.json` and prepares the files
//! their vsearch command lines read (`rscompare fetch`).
//!
//! Each dataset goes to `DIR/<id>/` (DIR defaults to `$RSCOMPARE_DATA`,
//! else `$TMPDIR/rsearch-data`, where the `datasets` source looks):
//!
//! - `DIR/<id>/download/`: the archive's files as published, checked
//!   against the manifest's size and MD5. A partial download is resumed; a
//!   file already complete is not fetched again.
//! - `DIR/<id>/`: the working files the command lines name. Each download
//!   is decompressed, and FASTQ is cut to its first N reads with
//!   `--max-reads` (pairs stay in step, since both mates are cut at N).
//!   Then come the dataset's `prepare` steps (pooled files, a simple exact
//!   dereplication, a UDB built by vsearch when `--vsearch` is given).
//!
//! Reference databases that a dataset's command lines use (their
//! `used_by`) are fetched with it. The working files are rebuilt only when
//! the downloads or `--max-reads` change. Downloads go through `curl`,
//! checksums through `md5` or `md5sum`, decompression through `gzip -dc`.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Component, Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread::sleep;
use std::time::Duration;

use serde_json::{Value, json};

const AGENT: &str = "rsearch-fetch-datasets/1";

/// What `rscompare fetch` was asked to do.
pub struct Options {
    pub ids: Vec<String>,
    pub dir: PathBuf,
    pub manifest: PathBuf,
    pub max_reads: Option<u64>,
    pub vsearch: Option<PathBuf>,
    pub list: bool,
    pub download_only: bool,
    pub no_references: bool,
    pub force: bool,
    pub retries: u32,
    pub quiet: bool,
}

/// `$RSCOMPARE_DATA`, else `$TMPDIR/rsearch-data`.
pub fn default_dir() -> PathBuf {
    match std::env::var_os("RSCOMPARE_DATA") {
        Some(d) if !d.is_empty() => PathBuf::from(d),
        _ => std::env::temp_dir().join("rsearch-data"),
    }
}

/// The manifest next to the crate.
pub fn default_manifest() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("datasets.json")
}

fn fail(msg: impl Into<String>) -> io::Error {
    io::Error::other(msg.into())
}

/// A size in bytes, as `12.3 MB`.
pub fn human(n: u64) -> String {
    if n < 1000 {
        return format!("{n} B");
    }
    let mut x = n as f64;
    for unit in ["kB", "MB", "GB"] {
        x /= 1000.0;
        if x < 1000.0 || unit == "GB" {
            return format!("{x:.1} {unit}");
        }
    }
    unreachable!()
}

/// The MD5 of a file, by `md5 -q` (macOS) or `md5sum`.
fn md5sum(path: &Path) -> io::Result<String> {
    for (prog, args) in [("md5", &["-q"][..]), ("md5sum", &[][..])] {
        if let Ok(out) = Command::new(prog)
            .args(args)
            .arg(path)
            .stderr(Stdio::null())
            .output()
            && out.status.success()
            && let Some(sum) = String::from_utf8_lossy(&out.stdout)
                .split_whitespace()
                .next()
        {
            return Ok(sum.to_lowercase());
        }
    }
    Err(fail("neither md5 nor md5sum works here"))
}

fn size_of(path: &Path) -> u64 {
    fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

/// Fetches `entry["url"]` into `dest`, resuming a `.part` file; checks size
/// and MD5. Returns whether it fetched something.
fn download(entry: &Value, dest: &Path, retries: u32, quiet: bool) -> io::Result<bool> {
    let size = entry["bytes"].as_u64();
    let md5 = entry["md5"].as_str().filter(|s| !s.is_empty());
    let url = entry["url"]
        .as_str()
        .ok_or_else(|| fail("a file without a url"))?;
    if dest.exists() {
        if let Some(size) = size
            && size_of(dest) != size
        {
            return Err(fail(format!(
                "{}: {} bytes, manifest says {size}; remove it to refetch",
                dest.display(),
                size_of(dest)
            )));
        }
        if let Some(md5) = md5
            && md5sum(dest)? != md5
        {
            return Err(fail(format!(
                "{}: MD5 differs from the manifest; remove it to refetch",
                dest.display()
            )));
        }
        return Ok(false);
    }
    let mut part = dest.as_os_str().to_owned();
    part.push(".part");
    let part = PathBuf::from(part);
    let mut attempt = 0;
    loop {
        let mut have = size_of(&part);
        if size.is_some_and(|s| have > s) {
            fs::remove_file(&part)?;
            have = 0;
        }
        if size == Some(have) {
            break;
        }
        // curl resumes from the part file's end (-C -), and starts over when
        // the server ignores the range
        let out = Command::new("curl")
            .args([
                "-fsSL",
                "--connect-timeout",
                "60",
                "--max-time",
                "0",
                "-A",
                AGENT,
                "-C",
                "-",
                "-o",
            ])
            .arg(&part)
            .arg(url)
            .stdout(Stdio::null())
            .stderr(Stdio::piped())
            .output()?;
        let error = if out.status.success() {
            if size.is_none() || size == Some(size_of(&part)) {
                break;
            }
            "connection closed early".to_string()
        } else {
            // 33: the server cannot resume; start over
            if out.status.code() == Some(33) {
                let _ = fs::remove_file(&part);
            }
            String::from_utf8_lossy(&out.stderr).trim().to_string()
        };
        if size_of(&part) > have {
            continue; // progress: resume at once
        }
        attempt += 1;
        if attempt > retries {
            return Err(fail(format!(
                "{url}: download failed after {retries} attempts without progress: {error}"
            )));
        }
        let wait = 60.min(1u64 << attempt);
        if !quiet {
            let name = dest.file_name().unwrap_or_default().to_string_lossy();
            eprintln!("  {name}: {error}; retry {attempt}/{retries} in {wait} s");
        }
        sleep(Duration::from_secs(wait));
    }
    let got = size_of(&part);
    if let Some(size) = size
        && got != size
    {
        return Err(fail(format!("{url}: got {got} bytes, expected {size}")));
    }
    if let Some(md5) = md5 {
        let actual = md5sum(&part)?;
        if actual != md5 {
            fs::remove_file(&part)?;
            return Err(fail(format!("{url}: MD5 {actual}, expected {md5}")));
        }
    }
    fs::rename(&part, dest)?;
    Ok(true)
}

/// A file's contents, through `gzip -dc` when it starts with gzip's magic
/// number.
fn open(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut magic = [0u8; 2];
    let n = File::open(path)?.read(&mut magic)?;
    if n == 2 && magic == [0x1f, 0x8b] {
        let child = Command::new("gzip")
            .arg("-dc")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()?;
        return Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            child.stdout.unwrap(),
        )));
    }
    Ok(Box::new(BufReader::with_capacity(
        1 << 20,
        File::open(path)?,
    )))
}

fn line(f: &mut dyn BufRead) -> io::Result<Vec<u8>> {
    let mut l = Vec::new();
    f.read_until(b'\n', &mut l)?;
    Ok(l)
}

fn trim_eol(mut l: &[u8]) -> &[u8] {
    while let [rest @ .., b'\r' | b'\n'] = l {
        l = rest;
    }
    l
}

fn trim(l: &[u8]) -> &[u8] {
    l.trim_ascii()
}

/// The next four-line FASTQ record, as its lines.
fn fastq_record(f: &mut dyn BufRead) -> io::Result<Option<[Vec<u8>; 4]>> {
    loop {
        let head = line(f)?;
        if head.is_empty() {
            return Ok(None);
        }
        if trim(&head).is_empty() {
            continue;
        }
        let rec = [head, line(f)?, line(f)?, line(f)?];
        if !rec[0].starts_with(b"@") || !rec[2].starts_with(b"+") {
            let shown = String::from_utf8_lossy(&rec[0][..rec[0].len().min(60)]).into_owned();
            return Err(fail(format!("not a four-line FASTQ record: {shown}")));
        }
        return Ok(Some(rec));
    }
}

/// The (label, sequence) pairs of a FASTA or FASTQ file.
fn sequences(path: &Path) -> io::Result<Vec<(Vec<u8>, Vec<u8>)>> {
    Ok(records(path)?
        .into_iter()
        .map(|(label, seq, _)| (label, seq))
        .collect())
}

/// A record: label, sequence and, from FASTQ, qualities.
type Rec = (Vec<u8>, Vec<u8>, Option<Vec<u8>>);

/// The records of a FASTA or FASTQ file.
fn records(path: &Path) -> io::Result<Vec<Rec>> {
    let mut f = open(path)?;
    let fastq = f.fill_buf()?.first() == Some(&b'@');
    let mut out = Vec::new();
    if fastq {
        while let Some(rec) = fastq_record(&mut *f)? {
            out.push((
                trim_eol(&rec[0][1..]).to_vec(),
                trim(&rec[1]).to_vec(),
                Some(trim(&rec[3]).to_vec()),
            ));
        }
        return Ok(out);
    }
    let mut cur: Option<Rec> = None;
    loop {
        let l = line(&mut *f)?;
        if l.is_empty() {
            break;
        }
        if l.starts_with(b">") {
            out.extend(cur.take());
            cur = Some((trim_eol(&l[1..]).to_vec(), Vec::new(), None));
        } else if let Some((_, s, _)) = &mut cur {
            s.extend_from_slice(trim(&l));
        }
    }
    out.extend(cur);
    Ok(out)
}

/// Writes through a `.tmp` file renamed into place at the end.
fn write_atomic(
    dst: &Path,
    body: impl FnOnce(&mut BufWriter<File>) -> io::Result<()>,
) -> io::Result<()> {
    let mut tmp = dst.as_os_str().to_owned();
    tmp.push(".tmp");
    let tmp = PathBuf::from(tmp);
    let mut out = BufWriter::with_capacity(1 << 20, File::create(&tmp)?);
    body(&mut out)?;
    out.flush()?;
    drop(out);
    fs::rename(&tmp, dst)
}

/// Decompresses `src` into `dst`, keeping the first `max_reads` FASTQ
/// records when given.
fn unpack(src: &Path, dst: &Path, max_reads: Option<u64>) -> io::Result<()> {
    let mut f = open(src)?;
    write_atomic(dst, |out| {
        match max_reads {
            Some(n) => {
                for _ in 0..n {
                    match fastq_record(&mut *f)? {
                        Some(rec) => rec.iter().try_for_each(|l| out.write_all(l))?,
                        None => break,
                    }
                }
            }
            None => {
                io::copy(&mut f, out)?;
            }
        }
        Ok(())
    })
}

/// A path of the manifest, `{dir}` the dataset's directory, relative to it.
fn resolve(path: &str, d: &Path) -> PathBuf {
    let p = PathBuf::from(path.replace("{dir}", &d.to_string_lossy()));
    let p = if p.is_absolute() { p } else { d.join(p) };
    // normalise as os.path.normpath does
    let mut out = PathBuf::new();
    for c in p.components() {
        match c {
            Component::ParentDir => {
                out.pop();
            }
            Component::CurDir => {}
            c => out.push(c),
        }
    }
    out
}

fn text<'a>(step: &'a Value, key: &str) -> io::Result<&'a str> {
    step[key]
        .as_str()
        .ok_or_else(|| fail(format!("a prepare step without {key}")))
}

fn op_cat(step: &Value, d: &Path) -> io::Result<()> {
    let inputs = step["inputs"]
        .as_array()
        .ok_or_else(|| fail("cat without inputs"))?;
    write_atomic(&resolve(text(step, "output")?, d), |out| {
        for name in inputs {
            let name = name
                .as_str()
                .ok_or_else(|| fail("cat input not a string"))?;
            io::copy(&mut File::open(resolve(name, d))?, out)?;
        }
        Ok(())
    })
}

fn op_fasta(step: &Value, d: &Path) -> io::Result<()> {
    let recs = sequences(&resolve(text(step, "input")?, d))?;
    write_atomic(&resolve(text(step, "output")?, d), |out| {
        for (label, seq) in recs {
            out.write_all(b">")?;
            out.write_all(&label)?;
            out.write_all(b"\n")?;
            out.write_all(&seq)?;
            out.write_all(b"\n")?;
        }
        Ok(())
    })
}

/// Exact dereplication, a stand-in for merging-free pipelines: upper case,
/// optionally truncated, filtered by length and Ns; sorted by decreasing
/// count, then first appearance; labels PREFIXn;size=N. With `"fastq":
/// true` it writes FASTQ, each sequence with the qualities of its first
/// read (from FASTQ input). It only makes realistic inputs; it is not
/// compared with anything.
fn op_derep(step: &Value, d: &Path) -> io::Result<()> {
    let fastq = step["fastq"].as_bool().unwrap_or(false);
    let trunc = step["truncate"]
        .as_u64()
        .filter(|&t| t > 0)
        .map(|t| t as usize);
    let lo = step["min_len"].as_u64().unwrap_or(1) as usize;
    let hi = step["max_len"]
        .as_u64()
        .filter(|&h| h > 0)
        .map(|h| h as usize);
    let max_n = step["max_n"].as_u64().map(|n| n as usize);
    let min_size = step["min_size"].as_u64().unwrap_or(1);
    let mut counts: HashMap<Vec<u8>, (u64, usize)> = HashMap::new();
    let mut quals: Vec<Vec<u8>> = Vec::new();
    for (_, seq, qual) in records(&resolve(text(step, "input")?, d))? {
        let mut seq = seq.to_ascii_uppercase();
        let mut qual = qual.unwrap_or_default();
        if let Some(t) = trunc {
            seq.truncate(t);
            qual.truncate(t);
        }
        if seq.len() < lo || hi.is_some_and(|h| seq.len() > h) {
            continue;
        }
        if max_n.is_some_and(|m| seq.iter().filter(|&&c| c == b'N').count() > m) {
            continue;
        }
        let order = counts.len();
        let entry = counts.entry(seq).or_insert((0, order));
        if entry.1 == quals.len() {
            quals.push(qual);
        }
        entry.0 += 1;
    }
    let mut keep: Vec<(Vec<u8>, (u64, usize))> = counts
        .into_iter()
        .filter(|(_, (c, _))| *c >= min_size)
        .collect();
    keep.sort_by_key(|(_, (c, o))| (std::cmp::Reverse(*c), *o));
    let prefix = step["prefix"].as_str().unwrap_or("U");
    write_atomic(&resolve(text(step, "output")?, d), |out| {
        for (k, (s, (c, o))) in keep.iter().enumerate() {
            if fastq {
                let q = &quals[*o];
                if q.len() != s.len() {
                    return Err(fail("derep to FASTQ needs FASTQ input"));
                }
                writeln!(out, "@{prefix}{};size={c}", k + 1)?;
                out.write_all(s)?;
                out.write_all(b"\n+\n")?;
                out.write_all(q)?;
            } else {
                writeln!(out, ">{prefix}{};size={c}", k + 1)?;
                out.write_all(s)?;
            }
            out.write_all(b"\n")?;
        }
        Ok(())
    })
}

/// Writes `benches.tsv`, the roles of the entry's `"benches"` field with
/// their resolved paths, which the benchmarks read
/// (`benches/common/data.rs`). Without the field, it removes the file.
/// Returns the roles whose files are missing.
fn write_bench_table(e: &Value, d: &Path) -> io::Result<Vec<String>> {
    let table = d.join("benches.tsv");
    let Some(roles) = e["benches"].as_object() else {
        if table.exists() {
            fs::remove_file(&table)?;
        }
        return Ok(Vec::new());
    };
    let mut missing = Vec::new();
    write_atomic(&table, |out| {
        for (role, path) in roles {
            let path = path
                .as_str()
                .ok_or_else(|| fail(format!("benches.{role} is not a path")))?;
            let path = resolve(path, d);
            if !path.exists() {
                missing.push(role.clone());
            }
            writeln!(out, "{role}\t{}", path.display())?;
        }
        Ok(())
    })?;
    Ok(missing)
}

/// Builds a UDB file with vsearch; false when it cannot (no vsearch, or
/// no input).
fn op_makeudb(step: &Value, d: &Path, vsearch: Option<&Path>) -> io::Result<bool> {
    let output = text(step, "output")?;
    let Some(vsearch) = vsearch else {
        eprintln!("  skipped {output}: needs --vsearch PATH (a UDB file is built by vsearch)");
        return Ok(false);
    };
    let src = resolve(text(step, "input")?, d);
    if !src.exists() {
        eprintln!("  skipped {output}: {} missing", src.display());
        return Ok(false);
    }
    let status = Command::new(vsearch)
        .arg("--makeudb_usearch")
        .arg(&src)
        .arg("--output")
        .arg(resolve(output, d))
        .args(["--threads", "1", "--quiet"])
        .status()?;
    if !status.success() {
        return Err(fail(format!(
            "vsearch --makeudb_usearch {} failed",
            src.display()
        )));
    }
    Ok(true)
}

/// The manifest's datasets.
pub fn load(path: &Path) -> io::Result<Vec<Value>> {
    let m: Value = serde_json::from_str(&fs::read_to_string(path)?).map_err(fail_json)?;
    let list = if m.is_object() {
        m["datasets"].clone()
    } else {
        m
    };
    list.as_array()
        .cloned()
        .ok_or_else(|| fail("the manifest has no datasets"))
}

fn fail_json(e: serde_json::Error) -> io::Error {
    fail(e.to_string())
}

fn id(e: &Value) -> &str {
    e["id"].as_str().unwrap_or("")
}

fn strs(v: &Value) -> Vec<&str> {
    v.as_array()
        .map(|a| a.iter().filter_map(Value::as_str).collect())
        .unwrap_or_default()
}

/// The chosen datasets and the references they use, references first (a
/// dataset's preparation may read them).
pub fn with_references(selected: &[&Value], every: &[Value]) -> Vec<Value> {
    let mut ids: Vec<&str> = selected.iter().map(|e| id(e)).collect();
    for e in every {
        if !ids.contains(&id(e)) && strs(&e["used_by"]).iter().any(|u| ids.contains(u)) {
            ids.push(id(e));
        }
    }
    let mut out: Vec<Value> = every
        .iter()
        .filter(|e| ids.contains(&id(e)))
        .cloned()
        .collect();
    out.sort_by_key(|e| e["kind"].as_str() != Some("reference"));
    out
}

fn files(e: &Value) -> &[Value] {
    e["files"].as_array().map(Vec::as_slice).unwrap_or(&[])
}

fn total(e: &Value) -> u64 {
    files(e)
        .iter()
        .map(|f| f["bytes"].as_u64().unwrap_or(0))
        .sum()
}

fn fetch_one(e: &Value, opts: &Options) -> io::Result<()> {
    let d = opts.dir.join(id(e));
    let dl = d.join("download");
    fs::create_dir_all(&dl)?;
    let fs_ = files(e);
    println!("{}: {} files, {}", id(e), fs_.len(), human(total(e)));
    for (k, f) in fs_.iter().enumerate() {
        let name = f["name"]
            .as_str()
            .ok_or_else(|| fail("a file without a name"))?;
        if download(f, &dl.join(name), opts.retries, opts.quiet)? && !opts.quiet {
            println!(
                "  [{}/{}] {name} {}",
                k + 1,
                fs_.len(),
                human(f["bytes"].as_u64().unwrap_or(0))
            );
        }
    }
    if opts.download_only {
        return Ok(());
    }
    let reads = if e["kind"].as_str() == Some("reference") {
        None
    } else {
        opts.max_reads
    };
    let names: serde_json::Map<String, Value> = fs_
        .iter()
        .map(|f| {
            (
                f["name"].as_str().unwrap_or("").to_string(),
                f["md5"].clone(),
            )
        })
        .collect();
    let prepare = e.get("prepare").cloned().unwrap_or(json!([]));
    let mut state = json!({
        "max_reads": reads,
        "files": names,
        "prepare": prepare,
        "vsearch": opts.vsearch.is_some(),
    });
    let stamp = d.join(".prepared.json");
    if stamp.exists() && !opts.force {
        let old: Option<Value> = fs::read_to_string(&stamp)
            .ok()
            .and_then(|t| serde_json::from_str(&t).ok());
        if old.as_ref() == Some(&state) {
            println!("  working files up to date");
            return bench_table(e, &d);
        }
    }
    for f in fs_ {
        let name = f["name"].as_str().unwrap_or("");
        let fastq = f["format"].as_str().unwrap_or("fastq") == "fastq";
        let unpacked = name.strip_suffix(".gz").unwrap_or(name);
        unpack(
            &dl.join(name),
            &d.join(unpacked),
            if fastq { reads } else { None },
        )?;
    }
    let mut complete = true;
    for step in prepare.as_array().map(Vec::as_slice).unwrap_or(&[]) {
        match step["op"].as_str() {
            Some("cat") => op_cat(step, &d)?,
            Some("fasta") => op_fasta(step, &d)?,
            Some("derep") => op_derep(step, &d)?,
            Some("makeudb") => complete &= op_makeudb(step, &d, opts.vsearch.as_deref())?,
            other => return Err(fail(format!("unknown prepare step {other:?}"))),
        }
    }
    if !complete {
        state["vsearch"] = Value::Null; // retry the skipped steps next time
    }
    fs::write(&stamp, state.to_string())?;
    let cut = reads
        .map(|r| format!(" (first {r} reads of each FASTQ)"))
        .unwrap_or_default();
    println!("  working files in {}{cut}", d.display());
    bench_table(e, &d)
}

/// Writes the benchmarks' table, warning of missing files.
fn bench_table(e: &Value, d: &Path) -> io::Result<()> {
    let missing = write_bench_table(e, d)?;
    if !missing.is_empty() {
        eprintln!(
            "  benches.tsv: no file yet for {} (a skipped prepare step?)",
            missing.join(", ")
        );
    }
    Ok(())
}

/// Lists, or downloads and prepares, the datasets.
pub fn run(opts: &Options) -> io::Result<()> {
    let every = load(&opts.manifest)?;
    if opts.list {
        for e in &every {
            let what = e["platform"].as_str().or(e["kind"].as_str()).unwrap_or("");
            println!(
                "{:28} {:3} files {:>9}  {what}",
                id(e),
                files(e).len(),
                human(total(e))
            );
        }
        let n: usize = every.iter().map(|e| files(e).len()).sum();
        let bytes: u64 = every.iter().map(total).sum();
        println!("{:28} {n:3} files {:>9}", "total", human(bytes));
        return Ok(());
    }
    let unknown: Vec<&str> = opts
        .ids
        .iter()
        .map(String::as_str)
        .filter(|i| !every.iter().any(|e| id(e) == *i))
        .collect();
    if !unknown.is_empty() {
        return Err(fail(format!(
            "unknown dataset {} (see --list)",
            unknown.join(", ")
        )));
    }
    let selected: Vec<&Value> = every
        .iter()
        .filter(|e| opts.ids.is_empty() || opts.ids.iter().any(|i| i == id(e)))
        .collect();
    let chosen = if opts.no_references {
        selected.into_iter().cloned().collect()
    } else {
        with_references(&selected, &every)
    };
    fs::create_dir_all(&opts.dir)?;
    for e in &chosen {
        fetch_one(e, opts)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes() {
        assert_eq!(human(999), "999 B");
        assert_eq!(human(1500), "1.5 kB");
        assert_eq!(human(238_708_802), "238.7 MB");
        assert_eq!(human(5_000_000_000_000), "5000.0 GB");
    }

    #[test]
    fn references_first() {
        let every = vec![
            json!({"id": "a"}),
            json!({"id": "r", "kind": "reference", "used_by": ["a"]}),
            json!({"id": "s", "kind": "reference", "used_by": ["b"]}),
        ];
        let got = with_references(&[&every[0]], &every);
        assert_eq!(got.iter().map(id).collect::<Vec<_>>(), ["r", "a"]);
    }

    #[test]
    fn prepare_steps() {
        let d = std::env::temp_dir().join(format!("rscompare-fetch-{}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        fs::write(
            d.join("a.fq"),
            "@r1 x\nacgt\n+\nIIII\n\n@r2\nACGTN\n+\nIIIII\n@r3\nACGT\n+\nIIII\n",
        )
        .unwrap();
        unpack(&d.join("a.fq"), &d.join("b.fq"), Some(2)).unwrap();
        assert_eq!(
            fs::read_to_string(d.join("b.fq")).unwrap(),
            "@r1 x\nacgt\n+\nIIII\n@r2\nACGTN\n+\nIIIII\n"
        );
        op_fasta(&json!({"input": "b.fq", "output": "b.fa"}), &d).unwrap();
        assert_eq!(
            fs::read_to_string(d.join("b.fa")).unwrap(),
            ">r1 x\nacgt\n>r2\nACGTN\n"
        );
        op_cat(
            &json!({"inputs": ["b.fa", "{dir}/b.fa"], "output": "c.fa"}),
            &d,
        )
        .unwrap();
        op_derep(
            &json!({"input": "c.fa", "output": "u.fa", "max_n": 0, "prefix": "X"}),
            &d,
        )
        .unwrap();
        assert_eq!(
            fs::read_to_string(d.join("u.fa")).unwrap(),
            ">X1;size=2\nACGT\n"
        );
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn derep_to_fastq_keeps_the_first_qualities() {
        let d = std::env::temp_dir().join(format!("rscompare-fetch-fq-{}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        fs::write(
            d.join("a.fq"),
            "@r1\nacgtA\n+\nABCDE\n@r2\nTT\n+\nII\n@r3\nACGTT\n+\nFGHIJ\n",
        )
        .unwrap();
        op_derep(
            &json!({"input": "a.fq", "output": "u.fq", "truncate": 4, "fastq": true}),
            &d,
        )
        .unwrap();
        assert_eq!(
            fs::read_to_string(d.join("u.fq")).unwrap(),
            "@U1;size=2\nACGT\n+\nABCD\n@U2;size=1\nTT\n+\nII\n"
        );
        fs::remove_dir_all(&d).unwrap();
    }

    #[test]
    fn bench_table() {
        let d = std::env::temp_dir().join(format!("rscompare-fetch-tsv-{}", std::process::id()));
        let data = d.join("set");
        fs::create_dir_all(&data).unwrap();
        fs::write(data.join("r.fq"), "").unwrap();
        fs::create_dir_all(d.join("ref")).unwrap();
        fs::write(d.join("ref/x.fa"), "").unwrap();
        let entry = json!({"benches": {
            "reads": "r.fq",
            "db": "{dir}/../ref/x.fa",
            "mate": "absent.fq"
        }});
        let missing = write_bench_table(&entry, &data).unwrap();
        assert_eq!(missing, ["mate"]);
        let table = fs::read_to_string(data.join("benches.tsv")).unwrap();
        assert_eq!(
            table,
            format!(
                "db\t{}\nmate\t{}\nreads\t{}\n",
                d.join("ref/x.fa").display(),
                data.join("absent.fq").display(),
                data.join("r.fq").display()
            )
        );
        // an entry without benches leaves no table behind
        write_bench_table(&json!({}), &data).unwrap();
        assert!(!data.join("benches.tsv").exists());
        fs::remove_dir_all(&d).unwrap();
    }
}
