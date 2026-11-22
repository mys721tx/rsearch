//! Real datasets as benchmark input, shared by the benchmarks of both
//! crates through each `benches/<tool>/input.rs`
//! (`#[path = "../common/data.rs"] pub mod data;`). Only std.
//!
//! With `RSEARCH_BENCH_DATA` unset, nothing here changes a benchmark. It
//! writes and times its synthetic input as before. Set to a dataset id of
//! `compare/datasets.json`, the benchmarks read that dataset's prepared
//! files instead, as `rscompare fetch <id>` leaves them in
//! `$RSCOMPARE_DATA/<id>/` (default `$TMPDIR/rsearch-data/<id>/`):
//!
//! - `benches.tsv` there, written by `rscompare fetch` from the entry's
//!   `"benches"` field, maps roles to files: `reads` (FASTQ, Phred+33),
//!   `mate` (the second mates of `reads`), `uniques` (dereplicated FASTA
//!   sorted by `;size=`), `uniques_fq` (the same as FASTQ), `long` (long
//!   dereplicated sequences), `db` (a reference with `;tax=`), `sff` and
//!   `phred64`. A key `<bench>.<role>` overrides a role for one bench.
//! - Each bench derives what else it needs from those files, as it does
//!   from its synthetic reads (the FASTA of them, Phred+64, label lists).
//! - Criterion groups get `@<id>` (`derep/binary@geo-miseq-16s-v3v4`). Real
//!   timings thus never mix with the synthetic history. Throughputs count
//!   the records of the real files.
//! - A bench whose roles the dataset lacks (pairs of a single-end
//!   dataset) says so and times nothing.

#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// The environment variable naming the dataset.
pub const VARIABLE: &str = "RSEARCH_BENCH_DATA";
/// The table of roles in a dataset's directory.
pub const TABLE: &str = "benches.tsv";

/// The data directory of `rscompare fetch`.
pub fn root() -> PathBuf {
    if let Ok(d) = std::env::var("RSCOMPARE_DATA")
        && !d.is_empty()
    {
        return PathBuf::from(d);
    }
    let tmp = std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".into());
    Path::new(&tmp).join("rsearch-data")
}

/// The dataset id, if one is set.
pub fn id() -> Option<String> {
    std::env::var(VARIABLE).ok().filter(|v| !v.is_empty())
}

/// A dataset's files by role.
#[derive(Debug)]
pub struct Dataset {
    pub id: String,
    pub dir: PathBuf,
    pub files: Vec<(String, PathBuf)>,
}

impl Dataset {
    /// The file of `role` for `bench`: `<bench>.<role>`, else `<role>`.
    pub fn lookup(&self, bench: &str, role: &str) -> Option<&Path> {
        let specific = format!("{bench}.{role}");
        let find = |key: &str| {
            self.files
                .iter()
                .find(|(k, _)| k == key)
                .map(|(_, p)| p.as_path())
        };
        find(&specific).or_else(|| find(role))
    }
}

/// Reads a table of `role<TAB>path` lines; paths are relative to `dir`.
pub fn parse_table(text: &str, dir: &Path) -> Vec<(String, PathBuf)> {
    text.lines()
        .map(str::trim_end)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter_map(|l| l.split_once('\t'))
        .map(|(k, p)| (k.to_string(), dir.join(p)))
        .collect()
}

/// The dataset named by `RSEARCH_BENCH_DATA`, or `None` for the synthetic
/// input. Panics when it is set but not fetched.
pub fn dataset() -> Option<Dataset> {
    let id = id()?;
    let dir = root().join(&id);
    let table = dir.join(TABLE);
    let text = std::fs::read_to_string(&table).unwrap_or_else(|e| {
        panic!(
            "{VARIABLE}={id}: cannot read {}: {e}; run `cargo run -p rscompare \
             --release -- fetch {id}` first",
            table.display()
        )
    });
    Some(Dataset {
        files: parse_table(&text, &dir),
        id,
        dir,
    })
}

/// The first of `roles` that `dataset` lacks for `bench`, if any.
pub fn missing_in(dataset: Option<&Dataset>, bench: &str, roles: &[&str]) -> Option<String> {
    let d = dataset?;
    roles
        .iter()
        .find(|r| d.lookup(bench, r).is_none())
        .map(|r| r.to_string())
}

/// Whether `bench` has to skip, because a dataset is set and lacks one of
/// `roles`. Says why on standard error.
pub fn skip(bench: &str, roles: &[&str]) -> bool {
    let d = dataset();
    match missing_in(d.as_ref(), bench, roles) {
        Some(role) => {
            let id = d.map(|d| d.id).unwrap_or_default();
            eprintln!("{bench}: dataset {id} has no `{role}` file; nothing to time");
            true
        }
        None => false,
    }
}

/// The file of `role` for `bench` in the dataset, or `None` without one.
/// Panics if the dataset lacks it (check with [`skip`] first) or the file
/// is not there.
pub fn file(bench: &str, role: &str) -> Option<PathBuf> {
    let d = dataset()?;
    let path = d
        .lookup(bench, role)
        .unwrap_or_else(|| panic!("{bench}: dataset {} has no `{role}` file", d.id))
        .to_path_buf();
    assert!(
        path.exists(),
        "{bench}: {} is missing; run `rscompare fetch {}`",
        path.display(),
        d.id
    );
    Some(path)
}

/// The file of an optional `role` for `bench`: `None` without a dataset,
/// or when the dataset has no such file.
pub fn optional(bench: &str, role: &str) -> Option<PathBuf> {
    optional_in(dataset().as_ref(), bench, role)
}

pub fn optional_in(dataset: Option<&Dataset>, bench: &str, role: &str) -> Option<PathBuf> {
    dataset?
        .lookup(bench, role)
        .filter(|p| p.exists())
        .map(Path::to_path_buf)
}

/// The text of `role` for `bench`, or `None` without a dataset.
pub fn text(bench: &str, role: &str) -> Option<String> {
    file(bench, role).map(|p| read(&p))
}

/// Reads a file of text. Panics on failure.
pub fn read(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The environment variable listing the thread counts of the threaded
/// tools' benchmarks, such as `1,4,10`.
pub const THREADS: &str = "RSEARCH_BENCH_THREADS";

/// The thread counts a threaded tool is timed at: those of
/// `RSEARCH_BENCH_THREADS`, else one and one per core.
pub fn thread_counts() -> Vec<usize> {
    if let Ok(list) = std::env::var(THREADS) {
        let counts: Vec<usize> = list
            .split(',')
            .map(|n| {
                n.trim()
                    .parse()
                    .unwrap_or_else(|_| panic!("{THREADS}: {list}"))
            })
            .collect();
        assert!(counts.iter().all(|&n| n > 0), "{THREADS}: {list}");
        return counts;
    }
    let cores = std::thread::available_parallelism().map_or(1, |n| n.get());
    if cores > 1 { vec![1, cores] } else { vec![1] }
}

/// A criterion group name: `name`, or `name@<id>` with a dataset.
pub fn group(name: &str) -> String {
    group_for(name, id().as_deref())
}

pub fn group_for(name: &str, id: Option<&str>) -> String {
    match id {
        Some(id) => format!("{name}@{id}"),
        None => name.to_string(),
    }
}

/// The throughput of a group: `synthetic` without a dataset, else the
/// records of `path`.
pub fn count(path: impl AsRef<Path>, synthetic: usize) -> u64 {
    if id().is_some() {
        records(&read(path.as_ref())) as u64
    } else {
        synthetic as u64
    }
}

/// A file in the target's temporary directory for a bench's derived
/// input, named by the dataset: `<bench>-<id>-<name>`.
pub fn tmp(bench: &str, name: &str) -> PathBuf {
    tmp_for(bench, name, &id().unwrap_or_default())
}

pub fn tmp_for(bench: &str, name: &str, id: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("{bench}-{id}-{name}"))
}

/// A file a bench writes (hits, tables): `path` itself without a dataset,
/// else with `<id>-` before its name. Real runs then never overwrite the
/// synthetic runs' files.
pub fn output(path: PathBuf) -> PathBuf {
    output_for(path, id().as_deref())
}

pub fn output_for(path: PathBuf, id: Option<&str>) -> PathBuf {
    match (id, path.file_name()) {
        (Some(id), Some(name)) => path.with_file_name(format!("{id}-{}", name.to_string_lossy())),
        _ => path,
    }
}

/// Writes `text` to the temporary file `name` of `bench` and returns its
/// path.
pub fn write(bench: &str, name: &str, text: impl AsRef<[u8]>) -> PathBuf {
    let path = tmp(bench, name);
    std::fs::write(&path, text).unwrap();
    path
}

/// A FASTA or FASTQ record.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub label: String,
    pub seq: String,
    pub qual: Option<String>,
}

/// The records of FASTA (lines joined) or four-line FASTQ text.
pub fn parse(text: &str) -> Vec<Record> {
    let mut lines = text.lines().map(|l| l.trim_end_matches('\r'));
    let mut out = Vec::new();
    if text.starts_with('@') {
        while let Some(head) = lines.next() {
            if head.is_empty() {
                continue;
            }
            let seq = lines.next().unwrap_or_default();
            lines.next();
            let qual = lines.next().unwrap_or_default();
            out.push(Record {
                label: head[1..].to_string(),
                seq: seq.to_string(),
                qual: Some(qual.to_string()),
            });
        }
    } else {
        for line in lines {
            if let Some(label) = line.strip_prefix('>') {
                out.push(Record {
                    label: label.to_string(),
                    seq: String::new(),
                    qual: None,
                });
            } else if let Some(r) = out.last_mut() {
                r.seq.push_str(line.trim());
            }
        }
    }
    out
}

/// The number of records of FASTA or four-line FASTQ text.
pub fn records(text: &str) -> usize {
    if text.starts_with('@') {
        text.lines().filter(|l| !l.is_empty()).count() / 4
    } else {
        text.lines().filter(|l| l.starts_with('>')).count()
    }
}

/// Records as unwrapped FASTA.
pub fn to_fasta(records: &[Record]) -> String {
    records
        .iter()
        .map(|r| format!(">{}\n{}\n", r.label, r.seq))
        .collect()
}

/// Records as FASTQ, `+` lines bare. Panics on a record without
/// qualities.
pub fn to_fastq(records: &[Record]) -> String {
    records
        .iter()
        .map(|r| {
            let q = r.qual.as_deref().expect("a FASTA record written as FASTQ");
            format!("@{}\n{}\n+\n{q}\n", r.label, r.seq)
        })
        .collect()
}

/// Four-line Phred+33 FASTQ text re-encoded as Phred+64.
pub fn phred64(fastq: &str) -> String {
    fastq
        .lines()
        .enumerate()
        .map(|(i, line)| {
            if i % 4 == 3 {
                line.bytes().map(|q| char::from(q + 31)).collect::<String>() + "\n"
            } else {
                format!("{line}\n")
            }
        })
        .collect()
}

/// The distinct sequences of `records`, in first appearance.
pub fn distinct(records: &[Record]) -> Vec<String> {
    let mut seen = std::collections::HashSet::new();
    records
        .iter()
        .filter(|r| seen.insert(r.seq.as_str()))
        .map(|r| r.seq.clone())
        .collect()
}

/// A database for exact searches of `records`, as the generated ones are
/// made: their distinct sequences (`d<i>`), then every fourth of them again
/// (`dup<i>`). Some queries then hit several targets.
pub fn exact_db(records: &[Record]) -> String {
    let seqs = distinct(records);
    let mut db: String = seqs
        .iter()
        .enumerate()
        .map(|(i, s)| format!(">d{i}\n{s}\n"))
        .collect();
    for (i, s) in seqs.iter().enumerate().step_by(4) {
        db.push_str(&format!(">dup{i}\n{s}\n"));
    }
    db
}

/// The sum of the `;size=` annotations of FASTA or FASTQ text.
pub fn total_size(text: &str) -> u64 {
    parse(text)
        .iter()
        .filter_map(|r| {
            r.label
                .split(';')
                .find_map(|a| a.strip_prefix("size="))
                .and_then(|v| v.parse::<u64>().ok())
        })
        .sum()
}

/// Records with `;ee=` (from their Phred+33 qualities, four decimals) and
/// `;length=` appended to their labels.
pub fn with_ee_length(records: &[Record]) -> Vec<Record> {
    records
        .iter()
        .map(|r| {
            let ee: f64 = r
                .qual
                .as_deref()
                .unwrap_or_default()
                .bytes()
                .map(|q| 10f64.powf(-f64::from(q.saturating_sub(33)) / 10.0))
                .sum();
            let label = format!(
                "{};ee={ee:.4};length={}",
                r.label.trim_end_matches(';'),
                r.seq.len()
            );
            Record { label, ..r.clone() }
        })
        .collect()
}

/// Shuffles `items` with a seeded xorshift64*.
pub fn shuffle<T>(items: &mut [T], seed: u64) {
    let mut state = seed.max(1);
    let mut next = || {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        state.wrapping_mul(0x2545_f491_4f6c_dd1d)
    };
    for i in (1..items.len()).rev() {
        let j = (next() % (i as u64 + 1)) as usize;
        items.swap(i, j);
    }
}
