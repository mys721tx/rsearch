//! Checks that every tool's `--help` and its page in `docs/` agree.
//!
//! For every tool in `src/bin/`, `--help` succeeds and has a long
//! description that begins with the short one. Its long options are those
//! the page lists under "Options", no more and no fewer. Fixed defaults are
//! shown by clap rather than written in the help text. Every example of the
//! page, a `console` block of `$ ` commands and their output, gives that
//! output when run in a copy of `docs/examples/`. The examples of
//! `README.md` are run the same way.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The tools and their binaries, built for the tests.
const TOOLS: [(&str, &str); 40] = [
    ("allpr", env!("CARGO_BIN_EXE_allpr")),
    ("clmsa", env!("CARGO_BIN_EXE_clmsa")),
    ("clspl", env!("CARGO_BIN_EXE_clspl")),
    ("clust", env!("CARGO_BIN_EXE_clust")),
    ("conv", env!("CARGO_BIN_EXE_conv")),
    ("dbornt", env!("CARGO_BIN_EXE_dbornt")),
    ("dbsrch", env!("CARGO_BIN_EXE_dbsrch")),
    ("derep", env!("CARGO_BIN_EXE_derep")),
    ("eegrid", env!("CARGO_BIN_EXE_eegrid")),
    ("eestat", env!("CARGO_BIN_EXE_eestat")),
    ("exsrch", env!("CARGO_BIN_EXE_exsrch")),
    ("fqchar", env!("CARGO_BIN_EXE_fqchar")),
    ("fqstat", env!("CARGO_BIN_EXE_fqstat")),
    ("getseq", env!("CARGO_BIN_EXE_getseq")),
    ("hitfmt", env!("CARGO_BIN_EXE_hitfmt")),
    ("joinpa", env!("CARGO_BIN_EXE_joinpa")),
    ("lcatax", env!("CARGO_BIN_EXE_lcatax")),
    ("lcmask", env!("CARGO_BIN_EXE_lcmask")),
    ("mchime", env!("CARGO_BIN_EXE_mchime")),
    ("mergpa", env!("CARGO_BIN_EXE_mergpa")),
    ("mkudb", env!("CARGO_BIN_EXE_mkudb")),
    ("otutab", env!("CARGO_BIN_EXE_otutab")),
    ("redig", env!("CARGO_BIN_EXE_redig")),
    ("relab", env!("CARGO_BIN_EXE_relab")),
    ("rerep", env!("CARGO_BIN_EXE_rerep")),
    ("revcomp", env!("CARGO_BIN_EXE_revcomp")),
    ("scrmbl", env!("CARGO_BIN_EXE_scrmbl")),
    ("shufl", env!("CARGO_BIN_EXE_shufl")),
    ("sintax", env!("CARGO_BIN_EXE_sintax")),
    ("sortby", env!("CARGO_BIN_EXE_sortby")),
    ("sqfilt", env!("CARGO_BIN_EXE_sqfilt")),
    ("sqtrim", env!("CARGO_BIN_EXE_sqtrim")),
    ("subseq", env!("CARGO_BIN_EXE_subseq")),
    ("subsmp", env!("CARGO_BIN_EXE_subsmp")),
    ("syncpa", env!("CARGO_BIN_EXE_syncpa")),
    ("topn", env!("CARGO_BIN_EXE_topn")),
    ("uchime", env!("CARGO_BIN_EXE_uchime")),
    ("udbinfo", env!("CARGO_BIN_EXE_udbinfo")),
    ("udbstat", env!("CARGO_BIN_EXE_udbstat")),
    ("unoise", env!("CARGO_BIN_EXE_unoise")),
];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Runs a tool with one argument and returns its standard output, which
/// must have succeeded.
fn help(bin: &str, flag: &str) -> String {
    let output = Command::new(bin).arg(flag).output().unwrap();
    assert!(output.status.success(), "{bin} {flag} failed");
    String::from_utf8(output.stdout).unwrap()
}

/// The long options of a help text: those on the lines that introduce an
/// option, indented by two or six spaces, `--help` aside.
fn help_options(text: &str) -> BTreeSet<String> {
    let mut options = BTreeSet::new();
    for line in text.lines() {
        let indent = line.len() - line.trim_start().len();
        let line = line.trim_start();
        if indent > 6 || !line.starts_with('-') {
            continue;
        }
        let line = line.split_once(", ").map_or(line, |(_, rest)| rest);
        if let Some(rest) = line.strip_prefix("--") {
            let name: String = rest
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '-')
                .collect();
            if name != "help" {
                options.insert(name);
            }
        }
    }
    options
}

/// The long options of a page: those named in the head of each item of its
/// "## Options" section, a "- `" line and the lines it runs on to up to the
/// first "`: ", `--help` aside.
fn page_options(page: &str) -> BTreeSet<String> {
    let mut heads = Vec::new();
    let mut head: Option<String> = None;
    let mut in_options = false;
    for line in page.lines() {
        if line.starts_with("## ") {
            in_options = line == "## Options";
            continue;
        }
        if in_options && line.starts_with("- `") {
            head = Some(String::new());
        }
        if let Some(text) = head.as_mut() {
            text.push(' ');
            match format!("{line} ").split_once("`: ") {
                Some((before, _)) => {
                    text.push_str(before);
                    heads.push(head.take().unwrap());
                }
                None => text.push_str(line),
            }
        }
    }
    let mut options = BTreeSet::new();
    for head in heads {
        for piece in head.split("--").skip(1) {
            let name: String = piece
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '-')
                .collect();
            if !name.is_empty() && name != "help" {
                options.insert(name);
            }
        }
    }
    options
}

/// The examples of a page: the commands of its `console` blocks, each with
/// the output that follows it.
fn examples(page: &str) -> Vec<(String, String)> {
    let mut examples: Vec<(String, String)> = Vec::new();
    let mut in_console = false;
    for line in page.lines() {
        if !in_console {
            in_console = line == "```console";
            continue;
        }
        if line == "```" {
            in_console = false;
        } else if let Some(command) = line.strip_prefix("$ ") {
            examples.push((command.to_string(), String::new()));
        } else {
            let (_, output) = examples.last_mut().expect("a block starts with $");
            output.push_str(line);
            output.push('\n');
        }
    }
    examples
}

fn copy_dir(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        if entry.file_type().unwrap().is_file() {
            fs::copy(entry.path(), to.join(entry.file_name())).unwrap();
        }
    }
}

#[test]
fn every_tool_has_a_page() {
    let mut bins: Vec<String> = fs::read_dir(root().join("src/bin"))
        .unwrap()
        .map(|e| {
            let name = e.unwrap().file_name().into_string().unwrap();
            name.strip_suffix(".rs").unwrap().to_string()
        })
        .collect();
    bins.sort();
    let tools: Vec<String> = TOOLS.iter().map(|(t, _)| t.to_string()).collect();
    assert_eq!(bins, tools, "the tools tested are those of src/bin");
    for tool in &tools {
        assert!(
            root().join("docs").join(format!("{tool}.md")).is_file(),
            "docs/{tool}.md is missing"
        );
    }
}

#[test]
fn help_has_a_long_description_beginning_with_the_short_one() {
    for (tool, bin) in TOOLS {
        let long = help(bin, "--help");
        let short = help(bin, "-h");
        let about = short.lines().next().unwrap();
        assert!(
            about.ends_with('.'),
            "{tool}: about does not end with a period"
        );
        let (description, _) = long.split_once("\nUsage:").unwrap();
        let paragraphs: Vec<&str> = description.trim().split("\n\n").collect();
        assert!(paragraphs.len() >= 2, "{tool}: no long_about");
        assert_eq!(paragraphs[0], about, "{tool}: long_about begins otherwise");
    }
}

#[test]
fn help_shows_fixed_defaults_through_clap() {
    for (tool, bin) in TOOLS {
        for line in help(bin, "--help").lines() {
            let Some((_, value)) = line.trim().rsplit_once("default to ") else {
                continue;
            };
            assert!(
                value.contains(' ') || ["stdin", "stdout", "all"].contains(&value),
                "{tool}: a fixed default in the help text: {line}"
            );
        }
    }
}

#[test]
fn page_options_are_those_of_help() {
    for (tool, bin) in TOOLS {
        let path = root().join("docs").join(format!("{tool}.md"));
        let page = fs::read_to_string(&path).unwrap_or_default();
        let in_help = help_options(&help(bin, "--help"));
        let in_page = page_options(&page);
        let missing: Vec<_> = in_help.difference(&in_page).collect();
        let extra: Vec<_> = in_page.difference(&in_help).collect();
        assert!(
            missing.is_empty() && extra.is_empty(),
            "{tool}: options missing from its page {missing:?}, not in --help {extra:?}"
        );
    }
}

#[test]
fn page_examples_give_their_output() {
    let bin_dir = Path::new(TOOLS[0].1).parent().unwrap();
    let path = format!(
        "{}:{}",
        bin_dir.display(),
        std::env::var("PATH").unwrap_or_default()
    );
    let pages = TOOLS
        .iter()
        .map(|(tool, _)| {
            (
                tool.to_string(),
                root().join("docs").join(format!("{tool}.md")),
            )
        })
        .chain([("README".to_string(), root().join("README.md"))]);
    for (tool, page) in pages {
        let page = fs::read_to_string(page).unwrap_or_default();
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("docs-{tool}"));
        let _ = fs::remove_dir_all(&dir);
        copy_dir(&root().join("docs/examples"), &dir);
        for (command, expected) in examples(&page) {
            let output = Command::new("sh")
                .arg("-c")
                .arg(format!("exec 2>&1; {command}"))
                .current_dir(&dir)
                .env("PATH", &path)
                .output()
                .unwrap();
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                expected,
                "{tool}: the output of `{command}`"
            );
        }
        fs::remove_dir_all(&dir).unwrap();
    }
}
