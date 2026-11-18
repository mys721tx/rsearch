//! Output labels: header annotations (`;size=`, `;ee=`, `;length=`) and
//! relabelling, following vsearch's rules.

use md5::Md5;
use sha1::{Digest, Sha1};

use crate::cluster::normalize;

/// The `[start, end)` byte range of the first `name` annotation of a
/// header: one that starts the header or follows a `;`, with a value of
/// digits (and dots, if `dots`) that ends the header or is followed by a
/// `;`. These are vsearch's rules.
fn find_attribute(header: &[u8], name: &[u8], dots: bool) -> Option<(usize, usize)> {
    let mut offset = 0;
    while offset + name.len() <= header.len() {
        let hit = header[offset..]
            .windows(name.len())
            .position(|w| w == name)?
            + offset;
        offset = hit;
        if hit > 0 && header[hit - 1] != b';' {
            offset += name.len() + 1;
            continue;
        }
        let value = hit + name.len();
        let length = header[value..]
            .iter()
            .take_while(|c| c.is_ascii_digit() || (dots && **c == b'.'))
            .count();
        if length == 0 {
            offset += name.len() + 1;
            continue;
        }
        let end = value + length;
        if end < header.len() && header[end] != b';' {
            offset += name.len() + length + 2;
            continue;
        }
        return Some((hit, end));
    }
    None
}

/// returns the abundance annotated in a whole header (`;size=N`), or 1. As
/// in vsearch, the annotation must end the header or be followed by a `;`.
/// A description after it hides it.
///
/// # Example
/// ```rust
/// # use rsearch::label::abundance;
/// assert_eq!(abundance("read1;size=4;"), 4);
/// assert_eq!(abundance("read1;size=4 1:N:0"), 1);
/// assert_eq!(abundance("read1 1:N:0;size=4"), 4);
/// ```
pub fn abundance(header: &str) -> u64 {
    find_attribute(header.as_bytes(), b"size=", false)
        .and_then(|(start, end)| header[start + 5..end].parse().ok())
        .filter(|&size| size > 0)
        .unwrap_or(1)
}

/// Formats expected errors with about four significant digits, as vsearch.
fn format_ee(ee: f64) -> String {
    let decimals = match ee {
        x if x < 1e-9 => 13,
        x if x < 1e-8 => 12,
        x if x < 1e-7 => 11,
        x if x < 1e-6 => 10,
        x if x < 1e-5 => 9,
        x if x < 1e-4 => 8,
        x if x < 1e-3 => 7,
        x if x < 1e-2 => 6,
        x if x < 1e-1 => 5,
        _ => 4,
    };
    format!("{ee:.decimals$}")
}

/// Appends the header with the given annotations removed, as vsearch does.
/// Each annotation goes with the `;` before it. Returns whether what was
/// written ends with a `;`.
fn push_stripped(out: &mut String, header: &str, size: bool, ee: bool, length: bool) -> bool {
    let bytes = header.as_bytes();
    let mut found: Vec<(usize, usize)> = [
        (size, &b"size="[..], false),
        (ee, &b"ee="[..], true),
        (length, &b"length="[..], false),
    ]
    .into_iter()
    .filter(|&(wanted, _, _)| wanted)
    .filter_map(|(_, name, dots)| find_attribute(bytes, name, dots))
    .collect();
    found.sort_unstable();

    let mut last = None;
    let mut emit = |chunk: &str, out: &mut String| {
        if let Some(c) = chunk.chars().last() {
            out.push_str(chunk);
            last = Some(c);
        }
    };
    if found.is_empty() {
        emit(header, out);
    } else {
        let mut prev_end = 0;
        for &(start, end) in &found {
            if start > prev_end + 1 {
                emit(&header[prev_end..start - 1], out);
            }
            prev_end = end;
        }
        if header.len() > prev_end + 1 {
            emit(&header[prev_end..], out);
        }
    }
    last == Some(';')
}

/// How output labels are rewritten.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum Relabel {
    /// Keep the input label.
    #[default]
    None,
    /// PREFIX1, PREFIX2, ... by output order.
    Prefix(String),
    /// The sequence itself.
    Sequence,
    /// The SHA1 digest of the normalized sequence.
    Sha1,
    /// The MD5 digest of the normalized sequence.
    Md5,
}

/// Output label settings.
///
/// # Example
/// ```rust
/// # use rsearch::label::{Label, Relabel};
/// let label = Label {
///     sizeout: true,
///     eeout: true,
///     ..Label::default()
/// };
/// // the old annotations are replaced
/// assert_eq!(
///     label.format("read1;size=3;ee=0.5", b"ACGT", 2, 1, Some(0.125)),
///     "read1;size=2;ee=0.1250"
/// );
///
/// let label = Label {
///     relabel: Relabel::Prefix("r".into()),
///     keep: true,
///     ..Label::default()
/// };
/// assert_eq!(label.format("read1 desc", b"ACGT", 1, 7, None), "r7 read1 desc");
/// ```
#[derive(Clone, Debug, Default)]
pub struct Label {
    pub relabel: Relabel,
    /// Append the input label after a new one.
    pub keep: bool,
    /// Write the abundance (`;size=`), replacing the input's.
    pub sizeout: bool,
    /// Remove the input's abundance.
    pub xsize: bool,
    /// Write the expected errors (`;ee=`), replacing the input's.
    pub eeout: bool,
    /// Remove the input's expected errors.
    pub xee: bool,
    /// Write the sequence length (`;length=`), replacing the input's.
    pub lengthout: bool,
    /// Remove the input's length.
    pub xlength: bool,
    /// Add a sample annotation (`;sample=`).
    pub sample: Option<String>,
    /// Append this text to the label.
    pub suffix: Option<String>,
}

impl Label {
    /// whether labels are written unchanged.
    pub fn is_identity(&self) -> bool {
        self.relabel == Relabel::None
            && !(self.sizeout
                || self.xsize
                || self.eeout
                || self.xee
                || self.lengthout
                || self.xlength
                || self.sample.is_some()
                || self.suffix.is_some())
    }

    /// returns the output label of a sequence with input label `header`,
    /// abundance `size`, output position `ordinal` (from 1) and expected
    /// errors `ee` (`None` for FASTA).
    pub fn format(
        &self,
        header: &str,
        seq: &[u8],
        size: u64,
        ordinal: u64,
        ee: Option<f64>,
    ) -> String {
        let mut out = String::with_capacity(header.len() + 32);
        let mut trailing = false;
        let hex = |digest: &[u8], out: &mut String| {
            for byte in digest {
                out.push_str(&format!("{byte:02x}"));
            }
        };
        match &self.relabel {
            Relabel::Prefix(prefix) => {
                out.push_str(prefix);
                out.push_str(&ordinal.to_string());
            }
            Relabel::Sequence => out.push_str(&String::from_utf8_lossy(seq)),
            Relabel::Sha1 => hex(&Sha1::digest(normalize(seq)), &mut out),
            Relabel::Md5 => hex(&Md5::digest(normalize(seq)), &mut out),
            Relabel::None => {
                trailing = push_stripped(
                    &mut out,
                    header,
                    self.xsize || self.sizeout,
                    self.xee || (self.eeout && ee.is_some()),
                    self.xlength || self.lengthout,
                );
            }
        }
        if let Some(suffix) = &self.suffix {
            out.push_str(suffix);
            if !suffix.is_empty() {
                trailing = suffix.ends_with(';');
            }
        }
        let mut separator = |out: &mut String| {
            if trailing {
                trailing = false;
            } else {
                out.push(';');
            }
        };
        if let Some(sample) = &self.sample {
            separator(&mut out);
            out.push_str("sample=");
            out.push_str(sample);
        }
        if self.sizeout {
            separator(&mut out);
            out.push_str(&format!("size={size}"));
        }
        if let (true, Some(ee)) = (self.eeout, ee) {
            separator(&mut out);
            out.push_str("ee=");
            out.push_str(&format_ee(ee));
        }
        if self.lengthout {
            separator(&mut out);
            out.push_str(&format!("length={}", seq.len()));
        }
        if self.keep && self.relabel != Relabel::None {
            out.push(' ');
            out.push_str(header);
        }
        out
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_format_ee() {
        assert_eq!(format_ee(0.0), "0.0000000000000");
        assert_eq!(format_ee(0.00012346), "0.0001235");
        assert_eq!(format_ee(0.000012346), "0.00001235");
        assert_eq!(format_ee(0.05), "0.05000");
        assert_eq!(format_ee(12.5), "12.5000");
    }

    fn ee_label(header: &str, ee: f64) -> String {
        let label = Label {
            eeout: true,
            ..Label::default()
        };
        label.format(header, b"", 1, 1, Some(ee))
    }

    #[test]
    fn test_ee_annotation() {
        assert_eq!(ee_label("a", 1.0), "a;ee=1.0000");
        assert_eq!(ee_label("a;", 1.0), "a;ee=1.0000");
        assert_eq!(ee_label("ee=2;a", 1.0), ";a;ee=1.0000");
        assert_eq!(ee_label("a;ee=2;", 1.0), "a;ee=1.0000");
        // not annotations: inside a word, or with an invalid value
        assert_eq!(ee_label("free=2", 1.0), "free=2;ee=1.0000");
        assert_eq!(ee_label("a;ee=x", 1.0), "a;ee=x;ee=1.0000");
    }

    #[test]
    fn test_strip() {
        let strip = |header: &str| {
            let mut out = String::new();
            push_stripped(&mut out, header, true, true, false);
            out
        };
        assert_eq!(strip("a;size=3;b"), "a;b");
        assert_eq!(strip("size=3;b"), ";b");
        assert_eq!(strip("a;size=3;"), "a");
        assert_eq!(strip("a;ee=0.1;size=3;b"), "a;b");
        assert_eq!(strip("a;length=3"), "a;length=3");
    }

    #[test]
    fn test_relabel() {
        let label = |relabel| Label {
            relabel,
            sizeout: true,
            ..Label::default()
        };
        assert_eq!(
            label(Relabel::Sequence).format("a", b"ACGT", 2, 1, None),
            "ACGT;size=2"
        );
        // digests of the normalized sequence: case and U/T do not matter
        let sha1 = label(Relabel::Sha1).format("a", b"acgu", 1, 1, None);
        assert_eq!(sha1, label(Relabel::Sha1).format("b", b"ACGT", 1, 1, None));
        assert_eq!(sha1.len(), 40 + ";size=1".len());
        let md5 = label(Relabel::Md5).format("a", b"ACGT", 1, 1, None);
        assert_eq!(md5, "f1f8f4bf413b16ad135722aa4591043e;size=1");
    }

    #[test]
    fn test_sample_suffix_length() {
        let label = Label {
            sample: Some("S1".into()),
            suffix: Some(";".into()),
            lengthout: true,
            xsize: true,
            ..Label::default()
        };
        assert_eq!(
            label.format("a;size=3;length=9", b"ACG", 3, 1, None),
            "a;sample=S1;length=3"
        );
        assert!(Label::default().is_identity());
    }
}
