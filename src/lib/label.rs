//! Output labels: header annotations (`;size=`, `;ee=`, `;length=`) and
//! relabelling, following vsearch's rules.

use md5::Md5;
use sha1::{Digest, Sha1};

use crate::cluster::normalize_base;

/// The `[start, end)` byte range of the first `name` annotation of a
/// header: one that starts the header or follows a `;`, with a value of
/// digits (and dots, if `dots`) that ends the header or is followed by a
/// `;`. These are vsearch's rules.
fn find_attribute(header: &[u8], name: &[u8], dots: bool) -> Option<(usize, usize)> {
    let is_value = |c: u8| c.is_ascii_digit() || (dots && c == b'.');
    find_attribute_with(header, name, is_value, false)
}

/// The `[start, end)` byte range of the first `name` annotation of a
/// header, as vsearch's `header_find_attribute`: one that starts the
/// header or follows a `;`, whose value is of the bytes `is_value` takes
/// (empty only if `allow_empty`), and ends the header or is followed by a
/// `;`.
pub(crate) fn find_attribute_with(
    header: &[u8],
    name: &[u8],
    is_value: impl Fn(u8) -> bool,
    allow_empty: bool,
) -> Option<(usize, usize)> {
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
        let length = header[value..].iter().take_while(|&&c| is_value(c)).count();
        if length == 0 && !allow_empty {
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
    size_annotation(header)
        .filter(|&size| size > 0)
        .unwrap_or(1)
}

/// returns the abundance annotated in a whole header as written, under
/// the rules of `abundance`, or `None`. A zero is returned as is, for
/// tools that refuse it as vsearch's dereplication does.
///
/// # Example
/// ```rust
/// # use rsearch::label::size_annotation;
/// assert_eq!(size_annotation("read1;size=0"), Some(0));
/// assert_eq!(size_annotation("read1;size=x"), None);
/// ```
pub fn size_annotation(header: &str) -> Option<u64> {
    find_attribute(header.as_bytes(), b"size=", false)
        .and_then(|(start, end)| header[start + 5..end].parse().ok())
}

/// An abundance annotation that vsearch refuses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AbundanceError {
    /// `;size=0`.
    Zero,
    /// More than a signed 64-bit integer holds.
    Range,
}

impl std::fmt::Display for AbundanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AbundanceError::Zero => write!(f, "invalid (zero) abundance annotation"),
            AbundanceError::Range => write!(f, "invalid (range error) abundance annotation"),
        }
    }
}

impl std::error::Error for AbundanceError {}

/// returns the abundance annotated in a whole header under the rules of
/// `abundance`, or 1, refusing a zero and a value above `i64::MAX` as
/// vsearch does whenever it reads abundances.
///
/// # Example
/// ```rust
/// # use rsearch::label::{AbundanceError, checked_abundance};
/// assert_eq!(checked_abundance("a;size=4"), Ok(4));
/// assert_eq!(checked_abundance("a"), Ok(1));
/// assert_eq!(checked_abundance("a;size=0"), Err(AbundanceError::Zero));
/// assert_eq!(
///     checked_abundance("a;size=9223372036854775808"),
///     Err(AbundanceError::Range)
/// );
/// ```
pub fn checked_abundance(header: &str) -> Result<u64, AbundanceError> {
    let Some((start, end)) = find_attribute(header.as_bytes(), b"size=", false) else {
        return Ok(1);
    };
    // the value is digits only, and parsing fails only on overflow
    match header[start + 5..end].parse::<i64>() {
        Ok(0) => Err(AbundanceError::Zero),
        Ok(size) => Ok(size as u64),
        Err(_) => Err(AbundanceError::Range),
    }
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

/// Appends the bytes as lower-case hexadecimal.
fn push_hex(bytes: &[u8], out: &mut String) {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    out.reserve(2 * bytes.len());
    for &byte in bytes {
        out.push(char::from(DIGITS[usize::from(byte >> 4)]));
        out.push(char::from(DIGITS[usize::from(byte & 15)]));
    }
}

/// Appends the hexadecimal digest of the normalized sequence, normalizing
/// it chunk by chunk into a buffer on the stack rather than a copy.
fn push_digest<D: Digest>(seq: &[u8], out: &mut String) {
    let mut hasher = D::new();
    let mut buffer = [0u8; 256];
    for chunk in seq.chunks(buffer.len()) {
        for (normalized, &base) in buffer.iter_mut().zip(chunk) {
            *normalized = normalize_base(base);
        }
        hasher.update(&buffer[..chunk.len()]);
    }
    push_hex(&hasher.finalize(), out);
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
    /// Add the number of sequences of a cluster (`;seqs=`).
    pub cluster_size: Option<u64>,
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
                || self.cluster_size.is_some()
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
        match &self.relabel {
            Relabel::Prefix(prefix) => {
                out.push_str(prefix);
                out.push_str(&ordinal.to_string());
            }
            Relabel::Sequence => out.push_str(&String::from_utf8_lossy(seq)),
            Relabel::Sha1 => push_digest::<Sha1>(seq, &mut out),
            Relabel::Md5 => push_digest::<Md5>(seq, &mut out),
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
        if let Some(n) = self.cluster_size {
            separator(&mut out);
            out.push_str(&format!("seqs={n}"));
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
        // longer than the normalizing buffer, and hashed in several chunks
        let long: Vec<u8> = b"acgu".repeat(200);
        let expected: Vec<u8> = b"ACGT".repeat(200);
        let digest = |seq: &[u8]| label(Relabel::Md5).format("a", seq, 1, 1, None);
        assert_eq!(digest(&long), digest(&expected));
        let mut whole = String::new();
        push_hex(&Md5::digest(&expected), &mut whole);
        assert_eq!(&digest(&expected)[..32], whole);
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
