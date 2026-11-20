//! Runs vsearch and rsearch pipelines.

use std::collections::HashSet;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::process::{Child, Command, Stdio};
use std::thread;
use std::time::{Duration, Instant};

use crate::translate::{Step, shell_quote};

/// How long one program may run.
pub const TIMEOUT: Duration = Duration::from_secs(1800);

/// The result of running vsearch or a pipeline.
#[derive(Debug, Clone, Default)]
pub struct Run {
    pub exit: i32,
    /// The failing step's (or vsearch's) standard error.
    pub stderr: String,
    /// Every step's standard error.
    pub all_stderr: String,
    pub seconds: f64,
    /// The first failing step, rendered.
    pub failed: Option<String>,
    /// Each step: its tool, exit code and time, or that it was skipped.
    pub steps: Vec<serde_json::Value>,
}

/// Waits for a child with a timeout, collecting its standard error; a
/// signal shows as 128 plus its number, as a shell reports it.
fn finish(mut child: Child) -> (i32, String) {
    let mut err = child.stderr.take();
    let reader = thread::spawn(move || {
        let mut buf = Vec::new();
        if let Some(e) = err.as_mut() {
            let _ = e.read_to_end(&mut buf);
        }
        buf
    });
    let start = Instant::now();
    let mut pause = Duration::from_micros(200);
    let status = loop {
        match child.try_wait() {
            Ok(Some(status)) => break Some(status),
            Ok(None) if start.elapsed() > TIMEOUT => {
                let _ = child.kill();
                let _ = child.wait();
                break None;
            }
            Ok(None) => {
                thread::sleep(pause);
                pause = (pause * 2).min(Duration::from_millis(20));
            }
            Err(_) => break None,
        }
    };
    let mut err = String::from_utf8_lossy(&reader.join().unwrap_or_default()).into_owned();
    let code = match status {
        None => {
            err.push_str("\nrscompare: timed out\n");
            137
        }
        Some(s) => s.code().unwrap_or_else(|| {
            use std::os::unix::process::ExitStatusExt;
            128 + s.signal().unwrap_or(0)
        }),
    };
    (code, err)
}

fn proc(argv: &[String], stdin: Option<&Path>, stdout: Option<&Path>, cwd: &Path) -> (i32, String) {
    let mut cmd = Command::new(&argv[0]);
    cmd.args(&argv[1..]).current_dir(cwd).stderr(Stdio::piped());
    match stdin.map(File::open) {
        None => cmd.stdin(Stdio::null()),
        Some(Ok(f)) => cmd.stdin(f),
        Some(Err(e)) => return (127, format!("rscompare: {e}\n")),
    };
    match stdout.map(File::create) {
        None => cmd.stdout(Stdio::null()),
        Some(Ok(f)) => cmd.stdout(f),
        Some(Err(e)) => return (127, format!("rscompare: {e}\n")),
    };
    match cmd.spawn() {
        Ok(child) => finish(child),
        Err(e) => (127, format!("rscompare: {}: {e}\n", argv[0])),
    }
}

/// Runs vsearch with `argv`.
pub fn run_vsearch(vsearch: &Path, argv: &[String], cwd: &Path) -> Run {
    let t = Instant::now();
    let mut all = vec![vsearch.to_string_lossy().into_owned()];
    all.extend(argv.iter().cloned());
    let (exit, err) = proc(&all, None, None, cwd);
    Run {
        exit,
        stderr: err.clone(),
        all_stderr: err,
        seconds: t.elapsed().as_secs_f64(),
        ..Run::default()
    }
}

/// Runs the steps in order. A failed step's output, and whatever is made
/// from it, is incomplete. A later step still reads it on its standard
/// input, as a shell pipeline would. A step that takes it as a file
/// argument (a database, a hit file, a second input) is skipped. The exit
/// code is the first failing step's. A failing program that is not
/// rsearch's (a decompressor) counts as 1, vsearch's for unreadable input.
pub fn run_steps(steps: &[Step], bindir: &Path, cwd: &Path, outdir: &Path) -> Run {
    let t0 = Instant::now();
    let prefix = format!("{}/", outdir.to_string_lossy().trim_end_matches('/'));
    let mut first: Option<(i32, String, &Step)> = None;
    let (mut errs, mut info) = (String::new(), Vec::new());
    let (mut made, mut tainted) = (HashSet::new(), HashSet::new());
    let abs = |f: &str| {
        if Path::new(f).is_absolute() {
            Path::new(f).to_path_buf()
        } else {
            cwd.join(f)
        }
    };
    for s in steps {
        let args_in: Vec<&String> = s.args.iter().filter(|a| made.contains(*a)).collect();
        let no_stdin = s.stdin.as_ref().is_some_and(|f| !abs(f).exists());
        if no_stdin || args_in.iter().any(|a| tainted.contains(*a)) {
            info.push(serde_json::json!({"tool": s.tool, "skipped": true}));
            let mut skipped: Vec<String> = s
                .args
                .iter()
                .filter(|a| a.starts_with(&prefix) && !made.contains(*a))
                .cloned()
                .collect();
            skipped.extend(s.stdout.iter().cloned());
            for f in skipped {
                made.insert(f.clone());
                tainted.insert(f);
            }
            continue;
        }
        let mut outs: Vec<String> = s
            .args
            .iter()
            .filter(|a| a.starts_with(&prefix) && !made.contains(*a))
            .cloned()
            .collect();
        outs.extend(s.stdout.iter().cloned());
        let t = Instant::now();
        let stdin = s.stdin.as_deref().map(abs);
        let stdout = s.stdout.as_deref().map(abs);
        let (mut code, err) = proc(&s.argv(bindir), stdin.as_deref(), stdout.as_deref(), cwd);
        if s.system && code != 0 {
            code = 1;
        }
        info.push(serde_json::json!({
            "tool": s.tool, "exit": code, "seconds": (t.elapsed().as_secs_f64() * 1e4).round() / 1e4,
        }));
        errs.push_str(&err);
        let reads_tainted = s.stdin.as_ref().is_some_and(|f| tainted.contains(f));
        for f in &outs {
            made.insert(f.clone());
        }
        if code != 0 || reads_tainted {
            tainted.extend(outs);
        }
        if code != 0 && first.is_none() {
            first = Some((code, err, s));
        }
    }
    let seconds = t0.elapsed().as_secs_f64();
    match first {
        None => Run {
            exit: 0,
            all_stderr: errs,
            seconds,
            steps: info,
            ..Run::default()
        },
        Some((code, err, s)) => Run {
            exit: code,
            stderr: err,
            all_stderr: errs,
            seconds,
            failed: Some(
                std::iter::once(s.tool.clone())
                    .chain(s.args.iter().cloned())
                    .collect::<Vec<_>>()
                    .join(" "),
            ),
            steps: info,
        },
    }
}

/// Runs the steps as a shell would, for timing. A step whose standard input
/// is the step before it's standard output (a temporary file used once)
/// reads it through a pipe. Any other step starts once the steps before it
/// have finished. Panics if a step fails.
pub fn run_piped(steps: &[Step], bindir: &Path) {
    let mut running: Vec<(String, Child)> = Vec::new();
    let wait = |running: &mut Vec<(String, Child)>| {
        for (cmd, mut child) in running.drain(..) {
            let status = child.wait().unwrap();
            assert!(status.success(), "{cmd} failed");
        }
    };
    let used_once = |file: &str| {
        steps
            .iter()
            .filter(|s| s.inputs().iter().any(|f| f == file))
            .count()
            == 1
    };
    for (k, s) in steps.iter().enumerate() {
        let argv = s.argv(bindir);
        let rendered = argv
            .iter()
            .map(|a| shell_quote(a))
            .collect::<Vec<_>>()
            .join(" ");
        let piped = k > 0
            && s.stdin.is_some()
            && s.stdin == steps[k - 1].stdout
            && used_once(s.stdin.as_deref().unwrap());
        let mut cmd = Command::new(&argv[0]);
        cmd.args(&argv[1..]).stderr(Stdio::null());
        if piped {
            let out = running.last_mut().and_then(|(_, c)| c.stdout.take());
            cmd.stdin(out.map_or(Stdio::null(), Stdio::from));
        } else {
            wait(&mut running);
            match &s.stdin {
                Some(f) => cmd.stdin(File::open(f).unwrap()),
                None => cmd.stdin(Stdio::null()),
            };
        }
        let next_piped = steps.get(k + 1).is_some_and(|n| {
            n.stdin.is_some() && n.stdin == s.stdout && used_once(n.stdin.as_deref().unwrap())
        });
        if next_piped {
            cmd.stdout(Stdio::piped());
        } else {
            match &s.stdout {
                Some(f) => cmd.stdout(File::create(f).unwrap()),
                None => cmd.stdout(Stdio::null()),
            };
        }
        running.push((rendered, cmd.spawn().unwrap()));
    }
    wait(&mut running);
}
