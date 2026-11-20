//! Builds vsearch from its source tree without autotools, for comparing
//! rsearch's output with vsearch's (`rscompare build-vsearch`).
//!
//! The architecture and the operating system pick vsearch's backends as
//! its `src/Makefile.am` does: `arch/x86_64` (SSE2 and SSSE3, one chosen at
//! run time), `arch/aarch64` (NEON), `arch/ppc64le` (AltiVec) or
//! `arch/generic`; `os/macos`, `os/freebsd`, or `os/linux` for Linux and
//! the other Unix systems, with `os/posix`. `HAVE_ZLIB_H` and
//! `HAVE_BZLIB_H` are defined when the compiler finds the headers, as
//! configure would, so that vsearch opens libz and libbz2 at run time to
//! read compressed input.

use std::collections::VecDeque;
use std::io;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{Arc, Mutex};

/// How to build.
#[derive(Debug, Clone)]
pub struct Options {
    /// The architecture, as `uname -m` or Rust names it (default: this one,
    /// or `VSEARCH_ARCH`).
    pub arch: String,
    /// The operating system, as `uname -s` or Rust names it (default: this
    /// one, or `VSEARCH_OS`).
    pub os: String,
    pub cxx: String,
    pub cc: String,
    pub jobs: usize,
    /// Print the chosen sources and flags, and build nothing.
    pub dry_run: bool,
}

impl Default for Options {
    fn default() -> Options {
        let env = |k: &str, d: &str| std::env::var(k).unwrap_or_else(|_| d.to_string());
        let arch = if std::env::consts::ARCH == "powerpc64" && cfg!(target_endian = "little") {
            "ppc64le"
        } else {
            std::env::consts::ARCH
        };
        Options {
            arch: env("VSEARCH_ARCH", arch),
            os: env("VSEARCH_OS", std::env::consts::OS),
            cxx: env("CXX", "c++"),
            cc: env("CC", "cc"),
            jobs: std::thread::available_parallelism().map_or(4, |n| n.get()),
            dry_run: false,
        }
    }
}

/// The platform's backends: the sources, those compiled with flags of
/// their own, the flags of every compile, and the libraries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Backends {
    pub sources: Vec<String>,
    pub simd: Vec<(String, String)>,
    pub flags: Vec<String>,
    pub libs: Vec<String>,
}

/// vsearch's backends for an architecture and an operating system.
pub fn backends(arch: &str, os: &str) -> io::Result<Backends> {
    let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
    let (mut sources, simd, flags) = match arch {
        "x86_64" | "amd64" => (
            s(&["arch/x86_64/cpu_features.cpp"]),
            vec![
                (
                    "arch/x86_64/SSE2/increment_counters.cpp".to_string(),
                    "-msse2".to_string(),
                ),
                (
                    "arch/x86_64/SSSE3/increment_counters.cpp".to_string(),
                    "-mssse3".to_string(),
                ),
            ],
            vec![],
        ),
        "arm64" | "aarch64" => (
            s(&[
                "arch/aarch64/increment_counters.cpp",
                "arch/generic/cpu_features.cpp",
            ]),
            vec![],
            vec![],
        ),
        "ppc64le" => (
            s(&[
                "arch/ppc64le/increment_counters.cpp",
                "arch/generic/cpu_features.cpp",
            ]),
            vec![],
            s(&["-mcpu=powerpc64le", "-maltivec"]),
        ),
        _ => (
            s(&[
                "arch/generic/increment_counters.cpp",
                "arch/generic/cpu_features.cpp",
            ]),
            vec![],
            vec![],
        ),
    };
    let posix = ["os/posix/dynlib_loader.cc", "os/posix/system.cc"];
    let mut libs = s(&["-lpthread"]);
    let lower = os.to_lowercase();
    if lower.starts_with("mingw")
        || lower.starts_with("msys")
        || lower.starts_with("cygwin")
        || lower.starts_with("windows")
    {
        return Err(io::Error::other("Windows builds need vsearch's own build"));
    }
    match lower.as_str() {
        "darwin" | "macos" => {
            sources.extend(s(&posix));
            sources.extend(s(&["os/macos/system_memory.cc", "os/macos/system_cpu.cc"]));
        }
        "freebsd" => {
            sources.extend(s(&posix));
            sources.extend(s(&[
                "os/freebsd/system_memory.cc",
                "os/freebsd/system_cpu.cc",
            ]));
        }
        _ => {
            // os/linux is vsearch's generic Unix backend
            sources.extend(s(&posix));
            sources.extend(s(&[
                "os/linux/system_memory.cc",
                "os/linux/system_cpu.cc",
                "os/linux/cgroup.cc",
            ]));
            if lower == "linux" {
                libs.push("-ldl".into());
            }
        }
    }
    Ok(Backends {
        sources,
        simd,
        flags,
        libs,
    })
}

/// The common sources of vsearch's `src/` tree: every C++ file but the
/// throwing fatal(), the library API, the platform backends and cityhash.
fn common_sources(src: &Path) -> io::Result<Vec<String>> {
    let mut out = Vec::new();
    let mut dirs = vec![src.to_path_buf()];
    while let Some(d) = dirs.pop() {
        for e in std::fs::read_dir(&d)? {
            let path = e?.path();
            if path.is_dir() {
                dirs.push(path);
                continue;
            }
            let rel = path
                .strip_prefix(src)
                .unwrap()
                .to_string_lossy()
                .replace('\\', "/");
            let cxx = rel.ends_with(".cc") || rel.ends_with(".cpp");
            let platform = rel.starts_with("arch/")
                || (rel.starts_with("os/") && rel.matches('/').count() >= 2);
            let left_out = rel.contains("fatal_throw")
                || rel.contains("vsearch_api")
                || rel == "vendored/city.cc";
            if cxx && !platform && !left_out {
                out.push(rel);
            }
        }
    }
    out.sort();
    Ok(out)
}

fn has_header(cxx: &str, header: &str) -> bool {
    let child = Command::new(cxx)
        .args(["-E", "-x", "c++", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn();
    let Ok(mut child) = child else { return false };
    if let Some(mut stdin) = child.stdin.take() {
        use std::io::Write;
        let _ = writeln!(stdin, "#include <{header}>");
    }
    child.wait().is_ok_and(|s| s.success())
}

fn obj(out: &Path, source: &str) -> PathBuf {
    out.join("obj")
        .join(format!("{}.o", source.replace('/', "_")))
}

fn compile(program: &str, args: &[String], cwd: &Path) -> io::Result<()> {
    let status = Command::new(program).args(args).current_dir(cwd).status()?;
    if status.success() {
        Ok(())
    } else {
        Err(io::Error::other(format!(
            "{program} {} failed",
            args.join(" ")
        )))
    }
}

/// Builds vsearch from the source tree `source` (the repository, holding
/// `src/`) into `out`, and returns the binary's path (`out/vsearch`).
pub fn build(source: &Path, out: &Path, opts: &Options) -> io::Result<PathBuf> {
    let src = source.join("src");
    let b = backends(&opts.arch, &opts.os)?;
    std::fs::create_dir_all(out.join("obj"))?;
    std::fs::create_dir_all(out.join("inc"))?;
    // what autotools would write: config.h for cityhash, and the package
    // name and version vsearch prints
    std::fs::write(out.join("inc/config.h"), "")?;
    std::fs::write(
        out.join("inc/package.h"),
        "#define PACKAGE \"vsearch\"\n#define PACKAGE_VERSION \"dev\"\n",
    )?;
    let mut flags: Vec<String> = ["-std=c++20", "-O2", "-DNDEBUG"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    flags.extend(b.flags.iter().cloned());
    flags.extend([
        "-include".to_string(),
        out.join("inc/package.h").to_string_lossy().into_owned(),
    ]);
    for (h, define) in [
        ("zlib.h", "-DHAVE_ZLIB_H=1"),
        ("bzlib.h", "-DHAVE_BZLIB_H=1"),
    ] {
        if has_header(&opts.cxx, h) {
            flags.push(define.into());
        } else {
            eprintln!("rscompare build-vsearch: no {h}, so no compressed input");
        }
    }
    let mut sources = common_sources(&src)?;
    sources.extend(b.sources.iter().cloned());
    if opts.dry_run {
        println!("arch={} os={}", opts.arch, opts.os);
        println!("flags: {}", flags.join(" "));
        println!("libs: {}", b.libs.join(" "));
        for s in sources
            .iter()
            .filter(|s| s.starts_with("arch/") || s.starts_with("os/"))
        {
            println!("{s}");
        }
        for (s, f) in &b.simd {
            println!("{s} ({f})");
        }
        return Ok(out.join("vsearch"));
    }
    for e in std::fs::read_dir(out.join("obj"))? {
        let p = e?.path();
        if p.extension().is_some_and(|x| x == "o") {
            std::fs::remove_file(p)?;
        }
    }
    // the common sources in parallel
    let queue = Arc::new(Mutex::new(sources.into_iter().collect::<VecDeque<_>>()));
    let failed = Arc::new(Mutex::new(None));
    std::thread::scope(|scope| {
        for _ in 0..opts.jobs.max(1) {
            let (queue, failed, flags) = (queue.clone(), failed.clone(), &flags);
            let (src, out, cxx) = (&src, out, &opts.cxx);
            scope.spawn(move || {
                loop {
                    let Some(s) = queue.lock().unwrap().pop_front() else {
                        break;
                    };
                    let mut args = flags.clone();
                    args.extend(
                        ["-fno-exceptions", "-I.", "-Ivendored", "-c"]
                            .iter()
                            .map(|x| x.to_string()),
                    );
                    args.extend([
                        s.clone(),
                        "-o".to_string(),
                        obj(out, &s).to_string_lossy().into_owned(),
                    ]);
                    if let Err(e) = compile(cxx, &args, src) {
                        *failed.lock().unwrap() = Some(e);
                        break;
                    }
                }
            });
        }
    });
    if let Some(e) = failed.lock().unwrap().take() {
        return Err(e);
    }
    for (s, f) in &b.simd {
        let mut args = flags.clone();
        args.extend([
            f.clone(),
            "-fno-exceptions".into(),
            "-I.".into(),
            "-c".into(),
            s.clone(),
            "-o".into(),
        ]);
        args.push(obj(out, s).to_string_lossy().into_owned());
        compile(&opts.cxx, &args, &src)?;
    }
    let inc = format!("-I{}", out.join("inc").to_string_lossy());
    let mut args = flags.clone();
    args.extend([
        inc,
        "-I.".into(),
        "-Ivendored".into(),
        "-c".into(),
        "vendored/city.cc".into(),
        "-o".into(),
    ]);
    args.push(out.join("obj/city.o").to_string_lossy().into_owned());
    compile(&opts.cxx, &args, &src)?;
    for c in ["md5", "sha1"] {
        let args: Vec<String> = [
            "-O2",
            "-c",
            &format!("vendored/{c}.c"),
            "-o",
            &out.join(format!("obj/{c}.o")).to_string_lossy(),
        ]
        .iter()
        .map(|x| x.to_string())
        .collect();
        compile(&opts.cc, &args, &src)?;
    }
    let mut objs: Vec<String> = std::fs::read_dir(out.join("obj"))?
        .filter_map(Result::ok)
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "o"))
        .map(|p| p.to_string_lossy().into_owned())
        .collect();
    objs.sort();
    let bin = out.join("vsearch");
    objs.extend(b.libs.iter().cloned());
    objs.extend(["-o".to_string(), bin.to_string_lossy().into_owned()]);
    compile(&opts.cxx, &objs, out)?;
    Ok(bin)
}

#[cfg(test)]
mod tests {
    use super::backends;

    fn chosen(arch: &str, os: &str) -> (Vec<String>, Vec<String>, Vec<String>) {
        let b = backends(arch, os).unwrap();
        let mut sources = b.sources;
        sources.extend(b.simd.iter().map(|(s, f)| format!("{s} {f}")));
        (sources, b.flags, b.libs)
    }

    #[test]
    fn as_makefile_am() {
        let (s, f, l) = chosen("x86_64", "Linux");
        assert!(s.contains(&"arch/x86_64/cpu_features.cpp".to_string()));
        assert!(s.contains(&"arch/x86_64/SSE2/increment_counters.cpp -msse2".to_string()));
        assert!(s.contains(&"arch/x86_64/SSSE3/increment_counters.cpp -mssse3".to_string()));
        assert!(s.contains(&"os/linux/cgroup.cc".to_string()));
        assert!(f.is_empty());
        assert_eq!(l, ["-lpthread", "-ldl"]);
        let (s, _, l) = chosen("arm64", "Darwin");
        assert!(s.contains(&"arch/aarch64/increment_counters.cpp".to_string()));
        assert!(s.contains(&"arch/generic/cpu_features.cpp".to_string()));
        assert!(s.contains(&"os/macos/system_cpu.cc".to_string()));
        assert_eq!(l, ["-lpthread"]);
        assert_eq!(chosen("aarch64", "macos").0, chosen("arm64", "Darwin").0);
        let (s, f, _) = chosen("ppc64le", "Linux");
        assert!(s.contains(&"arch/ppc64le/increment_counters.cpp".to_string()));
        assert_eq!(f, ["-mcpu=powerpc64le", "-maltivec"]);
        let (s, _, l) = chosen("amd64", "FreeBSD");
        assert!(s.contains(&"os/freebsd/system_memory.cc".to_string()));
        assert_eq!(l, ["-lpthread"]);
        let (s, _, l) = chosen("riscv64", "NetBSD");
        assert!(s.contains(&"arch/generic/increment_counters.cpp".to_string()));
        assert!(s.contains(&"os/linux/system_memory.cc".to_string()));
        assert_eq!(l, ["-lpthread"]);
        assert!(backends("x86_64", "MINGW64_NT").is_err());
    }
}
