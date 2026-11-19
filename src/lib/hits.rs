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
