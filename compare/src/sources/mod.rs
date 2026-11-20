//! Case sources. Each source writes its inputs with a generator
//! ([`crate::generators`]) from a seed and returns its cases: vsearch
//! command lines whose outputs go under [`OUT`](crate::cases::OUT).

use std::io;
use std::path::{Path, PathBuf};

use crate::cases::{Case, OUT};

pub mod analysis;
pub mod compressed;
pub mod datasets;
pub mod fuzz;
pub mod preprocess;

/// Where the programs are.
pub struct Env {
    /// rsearch's tools.
    pub bindir: PathBuf,
    pub vsearch: PathBuf,
}

/// Makes the cases of a seed in a directory.
pub type Make = fn(&Path, u64, &Env) -> io::Result<Vec<Case>>;

/// A source of cases: its default seeds, what it covers, how it makes them.
pub struct Source {
    pub name: &'static str,
    pub seeds: fn() -> Vec<u64>,
    pub about: &'static str,
    pub make: Make,
}

/// Every source, in the order `list` shows them.
pub fn all() -> Vec<Source> {
    let mut out = preprocess::sources();
    out.extend(analysis::sources());
    out.extend(fuzz::sources());
    out.extend(compressed::sources());
    out.extend(datasets::sources());
    out
}

pub fn find(name: &str) -> Option<Source> {
    all().into_iter().find(|s| s.name == name)
}

/// Seeds `a` to `b`.
pub fn seeds(a: u64, b: u64) -> Vec<u64> {
    (a..=b).collect()
}

/// Writes a generator's inputs.
pub fn generate(name: &str, dir: &Path, seed: u64, extra: &[&str]) -> io::Result<()> {
    crate::generators::write(name, dir, seed, extra)
}

/// Splits a command line as a shell does, quotes included, with `{O}`
/// standing for the output directory ([`OUT`]).
pub fn argv(text: &str) -> Vec<String> {
    shell_split(&text.replace("{O}", OUT))
}

/// Splits words as a POSIX shell does: single and double quotes, and
/// backslashes outside single quotes.
pub fn shell_split(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut word = String::new();
    let mut in_word = false;
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\'' => {
                in_word = true;
                for d in chars.by_ref() {
                    if d == '\'' {
                        break;
                    }
                    word.push(d);
                }
            }
            '"' => {
                in_word = true;
                while let Some(d) = chars.next() {
                    match d {
                        '"' => break,
                        '\\' => {
                            if let Some(e) = chars.next() {
                                word.push(e);
                            }
                        }
                        _ => word.push(d),
                    }
                }
            }
            '\\' => {
                in_word = true;
                if let Some(d) = chars.next() {
                    word.push(d);
                }
            }
            c if c.is_whitespace() => {
                if in_word {
                    out.push(std::mem::take(&mut word));
                    in_word = false;
                }
            }
            c => {
                in_word = true;
                word.push(c);
            }
        }
    }
    if in_word {
        out.push(word);
    }
    out
}

/// A file's text, or "".
pub fn read(path: &Path) -> String {
    std::fs::read_to_string(path).unwrap_or_default()
}

/// The words of a generator's `params` file.
pub fn params(dir: &Path) -> Vec<String> {
    read(&dir.join("params"))
        .split_whitespace()
        .map(str::to_string)
        .collect()
}

/// A path as text.
pub fn p(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

/// Strings from words.
pub fn strs(words: &[&str]) -> Vec<String> {
    words.iter().map(|w| w.to_string()).collect()
}

#[cfg(test)]
mod tests {
    use super::shell_split;

    #[test]
    fn splits() {
        assert_eq!(
            shell_split("a 'b c' \"d e\" f\\ g ''"),
            ["a", "b c", "d e", "f g", ""]
        );
        assert_eq!(
            shell_split("--label_suffix ';k=v;' x"),
            ["--label_suffix", ";k=v;", "x"]
        );
    }
}
