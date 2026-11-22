//! Timing a tool at each SIMD level this processor has (scalar, SSE2, and AVX2,
//! or scalar and NEON), chosen with `RSEARCH_SIMD`. The scalar level is plain
//! Rust that the compiler vectorises as it can. These show whether the
//! hand-written vector code still beats the compiler's.

#![allow(dead_code)]

use std::process::{Command, Stdio};

pub use rsearch::simd::Level;

/// the levels to time, those of `RSEARCH_BENCH_SIMD` (names such as
/// `scalar,avx2`), else every level this processor has.
pub fn levels() -> Vec<Level> {
    let available = rsearch::simd::available();
    match std::env::var("RSEARCH_BENCH_SIMD") {
        Ok(list) => list
            .split(',')
            .map(|name| {
                let level = Level::parse(name.trim())
                    .unwrap_or_else(|| panic!("RSEARCH_BENCH_SIMD: {list}"));
                assert!(available.contains(&level), "RSEARCH_BENCH_SIMD: {list}");
                level
            })
            .collect(),
        Err(_) => available,
    }
}

/// runs `program` with `args` at `level`, which must succeed.
pub fn run_at(level: Level, program: &str, args: &[&str]) {
    let status = Command::new(program)
        .env(rsearch::simd::VARIABLE, level.name())
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .unwrap();
    assert!(
        status.success(),
        "{program} {args:?} at {} failed",
        level.name()
    );
}
