//! Selecting records by their headers, as vsearch's `--fastx_getseq` and
//! `--fastx_getseqs` do.

use std::collections::HashSet;

/// Where a word may start and end: after and before a byte that is not an
/// ASCII letter or digit, or a `;` when the word is a field's value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Delimiter {
    NonAlphanumeric,
    Semicolon,
}

impl Delimiter {
    fn is(self, byte: u8) -> bool {
        match self {
            Delimiter::NonAlphanumeric => !byte.is_ascii_alphanumeric(),
            Delimiter::Semicolon => byte == b';',
        }
    }
}

/// whether `needle` occurs in `haystack` delimited at both ends, by a
/// delimiter or an end of the haystack. An empty needle is a word too, and
/// fits where two boundaries meet.
fn matches_delimited(haystack: &[u8], needle: &[u8], delimiter: Delimiter) -> bool {
    let n = needle.len();
    if n > haystack.len() {
        return false;
    }
    (0..=haystack.len() - n).any(|start| {
        let end = start + n;
        &haystack[start..end] == needle
            && (start == 0 || delimiter.is(haystack[start - 1]))
            && (end == haystack.len() || delimiter.is(haystack[end]))
    })
}

/// whether `needle` occurs in `haystack`, ignoring ASCII case. An empty
/// needle occurs in every haystack but an empty one, as with `std::search`.
fn contains_ignoring_case(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return !haystack.is_empty();
    }
    haystack
        .windows(needle.len())
        .any(|window| window.eq_ignore_ascii_case(needle))
}

/// What records are selected by.
#[derive(Clone, Debug)]
enum Rule {
    /// The whole header, ignoring ASCII case: upper-cased labels.
    Labels(HashSet<Vec<u8>>),
    /// A substring of the header, ignoring ASCII case.
    Substrings(Vec<Vec<u8>>),
    /// Delimited words: those without a delimiter, each of which can only
    /// be a whole token of the header, and the others.
    Words {
        tokens: HashSet<Vec<u8>>,
        delimited: Vec<Vec<u8>>,
        delimiter: Delimiter,
    },
}

/// Selects records by their headers.
///
/// # Example
/// ```rust
/// # use rsearch::select::Selector;
/// let label = Selector::labels(["seq1"], false);
/// assert!(label.matches("SEQ1"));
/// assert!(!label.matches("seq1 desc"));
///
/// let word = Selector::words(["s1"], Some("sample"));
/// assert!(word.matches("x;sample=s1;size=3"));
/// assert!(!word.matches("x;sample=s11"));
/// ```
#[derive(Clone, Debug)]
pub struct Selector {
    rule: Rule,
}

impl Selector {
    /// selects headers equal to one of `labels`, ignoring ASCII case, or
    /// that contain one of them with `substring`.
    pub fn labels<I, S>(labels: I, substring: bool) -> Selector
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let labels = labels.into_iter().map(|l| l.as_ref().to_vec());
        let rule = if substring {
            Rule::Substrings(labels.collect())
        } else {
            Rule::Labels(labels.map(|l| l.to_ascii_uppercase()).collect())
        };
        Selector { rule }
    }

    /// selects headers that contain one of `words` between delimiters, or
    /// with `field`, `field=word` between semicolons.
    pub fn words<I, S>(words: I, field: Option<&str>) -> Selector
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let (delimiter, prefix) = match field {
            Some(name) => (Delimiter::Semicolon, format!("{name}=").into_bytes()),
            None => (Delimiter::NonAlphanumeric, Vec::new()),
        };
        let (mut tokens, mut delimited) = (HashSet::new(), Vec::new());
        for word in words {
            let needle = [&prefix[..], word.as_ref()].concat();
            if needle.iter().any(|&b| delimiter.is(b)) {
                delimited.push(needle);
            } else {
                tokens.insert(needle);
            }
        }
        Selector {
            rule: Rule::Words {
                tokens,
                delimited,
                delimiter,
            },
        }
    }

    /// whether a header is selected.
    pub fn matches(&self, header: &str) -> bool {
        let header = header.as_bytes();
        match &self.rule {
            Rule::Labels(labels) => labels.contains(&header.to_ascii_uppercase()),
            Rule::Substrings(labels) => labels
                .iter()
                .any(|label| contains_ignoring_case(header, label)),
            Rule::Words {
                tokens,
                delimited,
                delimiter,
            } => {
                // a word without delimiters matches by being one whole token,
                // empty tokens included
                (!tokens.is_empty()
                    && header
                        .split(|&b| delimiter.is(b))
                        .any(|token| tokens.contains(token)))
                    || delimited
                        .iter()
                        .any(|needle| matches_delimited(header, needle, *delimiter))
            }
        }
    }
}

/// returns the labels of a labels file, one per line. Only a line feed ends
/// a line, and empty lines are skipped, as vsearch does.
///
/// # Example
/// ```rust
/// # use rsearch::select::labels_of;
/// assert_eq!(labels_of(b"a\n\nb\r\nc"), [&b"a"[..], b"b\r", b"c"]);
/// ```
pub fn labels_of(text: &[u8]) -> Vec<&[u8]> {
    text.split(|&b| b == b'\n')
        .filter(|line| !line.is_empty())
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_delimited() {
        let words = |w: &str| Selector::words([w], None);
        assert!(words("seq1").matches("y;tag=seq1"));
        assert!(!words("seq1").matches("seq10"));
        assert!(words("2").matches("s1-2.ok"));
        // an empty word sits where two boundaries meet
        assert!(words("").matches(""));
        assert!(words("").matches("a;;b"));
        assert!(!words("").matches("ab"));
        // a word with a delimiter is searched for
        assert!(words("1-2").matches("x 1-2 y"));
        assert!(!words("1-2").matches("s1-2.ok"));
    }

    #[test]
    fn test_substring() {
        let substr = |l: &str| Selector::labels([l], true);
        assert!(substr("EQ1").matches("seq10"));
        assert!(substr("").matches("a"));
        assert!(!substr("").matches(""));
    }
}
