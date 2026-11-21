//! What the benchmarks against vsearch share, through
//! `#[path = "../common/vsearch.rs"] mod vsearch;`.
//!
//! Each benchmark times vsearch commands (`vsearch/NAME`). For each, it
//! also times the rsearch pipeline the translator gives for that very
//! command (`pipeline/NAME`). The two thus always do the same job. They run
//! on the input files of rsearch's benchmark of the same tool. Both include
//! its generator (`benches/<tool>/input.rs`), and they use the same group
//! names.
//!
//! vsearch is `VSEARCH`, or else the build of `rscompare build-vsearch`;
//! without either a benchmark times nothing. rsearch's tools are those of
//! `cargo build --release`. Build them first:
//!
//! ```sh
//! cargo build --release
//! VSEARCH=/path/to/vsearch cargo bench -p rscompare --bench derep
//! ```

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use criterion::BenchmarkGroup;
use criterion::measurement::Measurement;

/// The target directory, the parent of `target/tmp`.
fn target() -> &'static Path {
    Path::new(env!("CARGO_TARGET_TMPDIR")).parent().unwrap()
}

/// The vsearch to time: `VSEARCH`, or else `compare/vsearch-build/vsearch`.
/// Prints why and returns `None` when there is neither.
pub fn vsearch() -> Option<String> {
    if let Ok(v) = std::env::var("VSEARCH") {
        return Some(v);
    }
    let built = Path::new(env!("CARGO_MANIFEST_DIR")).join("vsearch-build/vsearch");
    if built.exists() {
        return Some(built.to_string_lossy().into_owned());
    }
    eprintln!("no vsearch: set VSEARCH or run `rscompare build-vsearch`; nothing to time");
    None
}

/// The directory of rsearch's release tools. Panics if they are not built.
pub fn bindir() -> PathBuf {
    let dir = target().join("release");
    assert!(
        dir.join("conv").exists(),
        "rsearch's tools are not built: run `cargo build --release` first"
    );
    dir
}

/// The path of rsearch's tool `name`.
pub fn tool(name: &str) -> String {
    bindir().join(name).to_string_lossy().into_owned()
}

/// Runs `program` with `args`, discarding its output and errors. Panics if
/// it fails.
pub fn run<S: AsRef<str>>(program: &str, args: &[S]) {
    let status = Command::new(program)
        .args(args.iter().map(AsRef::as_ref))
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    let args: Vec<&str> = args.iter().map(AsRef::as_ref).collect();
    assert!(status.success(), "{program} {args:?} failed");
}

/// Registers `vsearch/NAME` (or `vsearch`), timing `vsearch` with `args`,
/// and `pipeline/NAME` (or `pipeline`), timing the rsearch pipeline they
/// translate to. Panics if they do not translate.
pub fn bench<M: Measurement, S: AsRef<str>>(
    group: &mut BenchmarkGroup<'_, M>,
    name: impl AsRef<str>,
    vsearch: &str,
    args: &[S],
) {
    let name = name.as_ref();
    let pipeline_name = name.replacen("vsearch", "pipeline", 1);
    let bindir = bindir();
    let dir = target()
        .join("tmp/rscompare-bench")
        .join(pipeline_name.replace('/', "-"));
    std::fs::create_dir_all(&dir).unwrap();
    let argv: Vec<&str> = args.iter().map(AsRef::as_ref).collect();
    let pipeline = rscompare::translate::translate(&argv, &dir.to_string_lossy())
        .unwrap_or_else(|e| panic!("vsearch {argv:?} does not translate: {e}"));
    group.bench_function(pipeline_name, |b| {
        b.iter(|| rscompare::execute::run_piped(&pipeline.steps, &bindir))
    });
    group.bench_function(name, |b| b.iter(|| run(vsearch, &argv)));
}
