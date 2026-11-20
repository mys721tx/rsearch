//! The rsearch pipeline a benchmark times against vsearch, shared by the
//! benchmarks through `#[path = "../common/translated.rs"] mod translated;`.
//!
//! For each vsearch command a benchmark times (`vsearch/NAME`), `bench`
//! times the pipeline compare/rscompare's translator gives for that very
//! command (`pipeline/NAME`), so the two always do the same job. The
//! single-tool benchmarks (`TOOL/NAME`) stay beside them for profiling.

use std::path::Path;

use criterion::BenchmarkGroup;
use criterion::measurement::Measurement;

/// Registers `pipeline/NAME` for the vsearch benchmark `vsearch/NAME` (or
/// `vsearch`, giving `pipeline`), timing the translation of its `vsearch`
/// arguments. Panics if they do not translate.
pub fn bench<M: Measurement, S: AsRef<str>>(
    group: &mut BenchmarkGroup<'_, M>,
    vsearch_name: impl AsRef<str>,
    vsearch: &[S],
) {
    let name = vsearch_name.as_ref().replacen("vsearch", "pipeline", 1);
    let bindir = Path::new(env!("CARGO_BIN_EXE_conv"))
        .parent()
        .unwrap()
        .to_path_buf();
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("rscompare-bench")
        .join(name.replace('/', "-"));
    std::fs::create_dir_all(&dir).unwrap();
    let argv: Vec<&str> = vsearch.iter().map(AsRef::as_ref).collect();
    let pipeline = rscompare::translate::translate(&argv, &dir.to_string_lossy())
        .unwrap_or_else(|e| panic!("vsearch {argv:?} does not translate: {e}"));
    group.bench_function(name, |b| {
        b.iter(|| rscompare::execute::run_piped(&pipeline.steps, &bindir))
    });
}
