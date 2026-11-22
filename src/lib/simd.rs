//! The vector instructions the tools run on, chosen once at run time.
//!
//! Every kernel that has vector versions asks [`level`] which one to run.
//! It is the best level the processor has, AVX2 or else SSE2 on x86-64 and
//! NEON on aarch64. `RSEARCH_SIMD` can choose a lower one, `scalar` on any
//! machine, so that each version can be run and compared with the others
//! on the same machine. Every level gives the same output, as every thread
//! count does. The binary is built for the baseline of its target. It runs
//! on any processor of that target.

use std::sync::OnceLock;

/// The environment variable that chooses a level below the best.
pub const VARIABLE: &str = "RSEARCH_SIMD";

/// A set of vector instructions.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Level {
    /// Plain Rust, on any machine.
    Scalar,
    /// 128-bit SSE2, part of every x86-64 processor.
    Sse2,
    /// 256-bit AVX2, on most x86-64 processors since 2013.
    Avx2,
    /// 128-bit NEON, part of every aarch64 processor.
    Neon,
}

/// Every level, by name.
const LEVELS: [Level; 4] = [Level::Scalar, Level::Sse2, Level::Avx2, Level::Neon];

impl Level {
    /// the name `RSEARCH_SIMD` takes.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::simd::Level;
    /// assert_eq!(Level::Avx2.name(), "avx2");
    /// assert_eq!(Level::parse("neon"), Some(Level::Neon));
    /// assert_eq!(Level::parse("sse4"), None);
    /// ```
    pub fn name(self) -> &'static str {
        match self {
            Level::Scalar => "scalar",
            Level::Sse2 => "sse2",
            Level::Avx2 => "avx2",
            Level::Neon => "neon",
        }
    }

    /// the level of a name.
    pub fn parse(name: &str) -> Option<Level> {
        LEVELS.into_iter().find(|level| level.name() == name)
    }
}

/// the levels this processor runs, the best last.
pub fn available() -> Vec<Level> {
    let mut levels = vec![Level::Scalar];
    #[cfg(target_arch = "x86_64")]
    {
        levels.push(Level::Sse2);
        if std::arch::is_x86_feature_detected!("avx2") {
            levels.push(Level::Avx2);
        }
    }
    #[cfg(target_arch = "aarch64")]
    levels.push(Level::Neon);
    levels
}

/// the level to run, the one `variable` names, which must be among
/// `available`, or else the best of them.
///
/// # Example
/// ```rust
/// # use rsearch::simd::{Level, choose};
/// let x86 = [Level::Scalar, Level::Sse2, Level::Avx2];
/// assert_eq!(choose(None, &x86), Ok(Level::Avx2));
/// assert_eq!(choose(Some("sse2"), &x86), Ok(Level::Sse2));
/// assert!(choose(Some("neon"), &x86).is_err());
/// assert!(choose(Some("fast"), &x86).is_err());
/// ```
pub fn choose(variable: Option<&str>, available: &[Level]) -> Result<Level, String> {
    let names = || {
        available
            .iter()
            .map(|level| level.name())
            .collect::<Vec<_>>()
            .join(", ")
    };
    let Some(name) = variable else {
        return Ok(*available.last().expect("scalar is always available"));
    };
    match Level::parse(name) {
        Some(level) if available.contains(&level) => Ok(level),
        Some(_) => Err(format!(
            "{VARIABLE}={name}: this processor does not have it (it has {})",
            names()
        )),
        None => Err(format!(
            "{VARIABLE}={name}: not a level (this processor has {})",
            names()
        )),
    }
}

/// the level the kernels run at, chosen on the first call. A value of
/// `RSEARCH_SIMD` that names no level, or one this processor lacks, stops
/// the tool with exit code 2, as an invalid argument does.
pub fn level() -> Level {
    static LEVEL: OnceLock<Level> = OnceLock::new();
    *LEVEL.get_or_init(|| {
        let variable = std::env::var(VARIABLE).ok();
        match choose(variable.as_deref(), &available()) {
            Ok(level) => level,
            Err(message) => {
                eprintln!("Error: {message}");
                std::process::exit(2);
            }
        }
    })
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn names_every_level_once() {
        for level in LEVELS {
            assert_eq!(Level::parse(level.name()), Some(level));
        }
        assert_eq!(Level::parse("AVX2"), None);
        assert_eq!(Level::parse(""), None);
    }

    #[test]
    fn has_scalar_and_its_target_baseline() {
        let levels = available();
        assert_eq!(levels[0], Level::Scalar);
        #[cfg(target_arch = "x86_64")]
        assert_eq!(levels[1], Level::Sse2);
        #[cfg(target_arch = "aarch64")]
        assert_eq!(levels, [Level::Scalar, Level::Neon]);
        assert!(levels.windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn chooses_the_best_or_the_one_named() {
        let neon = [Level::Scalar, Level::Neon];
        assert_eq!(choose(None, &neon), Ok(Level::Neon));
        assert_eq!(choose(Some("scalar"), &neon), Ok(Level::Scalar));
        let sse2 = [Level::Scalar, Level::Sse2];
        assert_eq!(
            choose(Some("avx2"), &sse2),
            Err("RSEARCH_SIMD=avx2: this processor does not have it (it has scalar, sse2)".into())
        );
        assert_eq!(
            choose(Some("Neon"), &neon),
            Err("RSEARCH_SIMD=Neon: not a level (this processor has scalar, neon)".into())
        );
    }

    #[test]
    fn runs_at_the_best_level_by_default() {
        if std::env::var(VARIABLE).is_err() {
            assert_eq!(level(), *available().last().unwrap());
        }
    }
}
