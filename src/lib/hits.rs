//! Hit records in vsearch's UC format, in which the tools pass hits on. It
//! holds both labels, the strand, the identity and the CIGAR, from which
//! hitfmt rebuilds every alignment. One line per hit, `H`, the target's
//! number (from 0), the query's length, the identity, the strand, two
//! zeros, the alignment as a CIGAR (or `=` for a perfect one), the query's
//! label and the target's; or `N` for a query without hits.

use std::io::{self, Write};

/// writes a hit.
///
/// # Example
/// ```rust
/// # use rsearch::hits::{write_hit, write_no_hit};
/// let mut out = Vec::new();
/// write_hit(&mut out, 3, 10, 100.0, false, b"=", b"q1", b"t3").unwrap();
/// write_no_hit(&mut out, b"q2").unwrap();
/// assert_eq!(
///     out,
///     b"H\t3\t10\t100.0\t+\t0\t0\t=\tq1\tt3\nN\t*\t*\t*\t.\t*\t*\t*\tq2\t*\n"
/// );
/// ```
#[allow(clippy::too_many_arguments)]
pub fn write_hit(
    out: &mut impl Write,
    target: u64,
    query_len: usize,
    identity: f64,
    minus: bool,
    cigar: &[u8],
    query: &[u8],
    target_label: &[u8],
) -> io::Result<()> {
    write!(
        out,
        "H\t{target}\t{query_len}\t{identity:.1}\t{}\t0\t0\t",
        if minus { '-' } else { '+' }
    )?;
    out.write_all(cigar)?;
    out.write_all(b"\t")?;
    out.write_all(query)?;
    out.write_all(b"\t")?;
    out.write_all(target_label)?;
    out.write_all(b"\n")
}

/// writes the record of a query without hits.
pub fn write_no_hit(out: &mut impl Write, query: &[u8]) -> io::Result<()> {
    out.write_all(b"N\t*\t*\t*\t.\t*\t*\t*\t")?;
    out.write_all(query)?;
    out.write_all(b"\t*\n")
}

/// A UC record: a hit, or a query without hits.
#[derive(Clone, Debug, PartialEq)]
pub struct UcRecord {
    /// `None` for an `N` record.
    pub hit: Option<UcHit>,
    /// The fields after the eighth: the query's label, a tab and the
    /// target's (`*` for an `N` record). Labels may hold tabs and are told
    /// apart with the target's label (see `split_labels`).
    pub labels: String,
}

/// The fields of an `H` record that hits are rebuilt from.
#[derive(Clone, Debug, PartialEq)]
pub struct UcHit {
    pub target: u64,
    pub query_len: usize,
    pub minus: bool,
    /// The CIGAR, or `=` for a perfect alignment.
    pub cigar: String,
}

impl UcRecord {
    /// parses a line of a UC file of hits, without its line feed.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::hits::UcRecord;
    /// let record = UcRecord::parse("H\t3\t10\t90.0\t-\t0\t0\t2I8M\tq 1\tt\t3").unwrap();
    /// let hit = record.hit.unwrap();
    /// assert_eq!((hit.target, hit.query_len, hit.minus), (3, 10, true));
    /// assert_eq!(record.labels, "q 1\tt\t3");
    /// assert!(UcRecord::parse("N\t*\t*\t*\t.\t*\t*\t*\tq\t*").unwrap().hit.is_none());
    /// assert!(UcRecord::parse("S\t0\t10\t*").is_err());
    /// ```
    pub fn parse(line: &str) -> Result<UcRecord, String> {
        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        let invalid = || format!("invalid hit record: {line}");
        if fields.len() < 9 {
            return Err(invalid());
        }
        let labels = fields[8].to_string();
        match fields[0] {
            "N" => Ok(UcRecord { hit: None, labels }),
            "H" => {
                let target = fields[1].parse().map_err(|_| invalid())?;
                let query_len = fields[2].parse().map_err(|_| invalid())?;
                let minus = match fields[4] {
                    "+" => false,
                    "-" => true,
                    _ => return Err(invalid()),
                };
                Ok(UcRecord {
                    hit: Some(UcHit {
                        target,
                        query_len,
                        minus,
                        cigar: fields[7].to_string(),
                    }),
                    labels,
                })
            }
            _ => Err(invalid()),
        }
    }

    /// the query's label, given the target's (`*` for an `N` record), or
    /// none if the record does not end with it.
    pub fn query_label(&self, target_label: &str) -> Option<&str> {
        self.labels
            .strip_suffix(target_label)
            .and_then(|rest| rest.strip_suffix('\t'))
    }
}
