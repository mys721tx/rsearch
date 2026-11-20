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

/// Reads a UC file of hits query by query: the consecutive `H` records of
/// one query label, or one `N` record. Labels are told apart with the
/// targets' labels, which the caller looks up by target number.
pub struct UcQueries<R> {
    lines: std::io::Lines<R>,
    pending: Option<UcRecord>,
}

impl<R: std::io::BufRead> UcQueries<R> {
    pub fn new(reader: R) -> UcQueries<R> {
        UcQueries {
            lines: reader.lines(),
            pending: None,
        }
    }

    fn record(&mut self) -> Option<Result<UcRecord, String>> {
        if let Some(record) = self.pending.take() {
            return Some(Ok(record));
        }
        let line = self.lines.next()?;
        Some(
            line.map_err(|e| format!("cannot read the hits: {e}"))
                .and_then(|line| UcRecord::parse(&line)),
        )
    }

    /// the next query's label and hits, none for an `N` record; `target`
    /// gives a target's label, or none beyond the database.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::hits::UcQueries;
    /// let uc = "H\t0\t4\t100.0\t+\t0\t0\t=\tq\ta\nH\t1\t4\t100.0\t+\t0\t0\t=\tq\tb\n\
    ///           N\t*\t*\t*\t.\t*\t*\t*\tr\t*\n";
    /// let labels = ["a", "b"];
    /// let target = |t: u64| labels.get(t as usize).copied();
    /// let mut queries = UcQueries::new(uc.as_bytes());
    /// let (label, hits) = queries.next_query(target).unwrap().unwrap();
    /// assert_eq!((label.as_str(), hits.len()), ("q", 2));
    /// let (label, hits) = queries.next_query(target).unwrap().unwrap();
    /// assert_eq!((label.as_str(), hits.len()), ("r", 0));
    /// assert!(queries.next_query(target).is_none());
    /// ```
    pub fn next_query<'a>(
        &mut self,
        target: impl Fn(u64) -> Option<&'a str>,
    ) -> Option<Result<(String, Vec<UcHit>), String>> {
        let label_of = |record: &UcRecord| -> Result<String, String> {
            let target_label = match &record.hit {
                None => "*",
                Some(hit) => target(hit.target)
                    .ok_or_else(|| format!("hit to target {}, beyond the database", hit.target))?,
            };
            record
                .query_label(target_label)
                .map(str::to_string)
                .ok_or_else(|| format!("labels do not match the database: {}", record.labels))
        };
        let first = match self.record()? {
            Ok(record) => record,
            Err(e) => return Some(Err(e)),
        };
        let label = match label_of(&first) {
            Ok(label) => label,
            Err(e) => return Some(Err(e)),
        };
        let mut hits = Vec::new();
        if let Some(hit) = first.hit {
            hits.push(hit);
            while let Some(next) = self.record() {
                let next = match next {
                    Ok(record) => record,
                    Err(e) => return Some(Err(e)),
                };
                let same = next.hit.is_some()
                    && match label_of(&next) {
                        Ok(l) => l == label,
                        Err(e) => return Some(Err(e)),
                    };
                if same {
                    hits.push(next.hit.unwrap());
                } else {
                    self.pending = Some(next);
                    break;
                }
            }
        }
        Some(Ok((label, hits)))
    }
}

/// An S or H record of clustering hits: the cluster of a record, its
/// label, and for a member its strand and its alignment to the centroid.
pub struct ClusterRecord<'a> {
    pub cluster: usize,
    pub label: &'a str,
    pub minus: bool,
    /// the CIGAR, `=` for a perfect alignment; empty for a centroid.
    pub cigar: &'a str,
}

impl<'a> ClusterRecord<'a> {
    /// parses an S or H record of the clusters so far, `clusters` of them:
    /// a centroid opens the next cluster, a member joins one of them. A C
    /// record gives `None`.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::hits::ClusterRecord;
    /// let s = ClusterRecord::parse("S\t0\t4\t*\t*\t*\t*\t*\ta\t*", 0).unwrap().unwrap();
    /// assert_eq!((s.cluster, s.label, s.cigar), (0, "a", ""));
    /// let h = ClusterRecord::parse("H\t0\t5\t80.0\t-\t0\t0\t4MD\tb\ta", 1).unwrap().unwrap();
    /// assert_eq!((h.cluster, h.label, h.minus, h.cigar), (0, "b", true, "4MD"));
    /// assert!(ClusterRecord::parse("C\t0\t2\t*\t*\t*\t*\t*\ta\t*", 1).unwrap().is_none());
    /// assert!(ClusterRecord::parse("H\t1\t5\t80.0\t+\t0\t0\t5M\tb\ta", 1).is_err());
    /// ```
    pub fn parse(line: &'a str, clusters: usize) -> Result<Option<ClusterRecord<'a>>, String> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 10 {
            return Err(format!("not a UC record: {line}"));
        }
        if fields[0] == "C" {
            return Ok(None);
        }
        let cluster: usize = fields[1]
            .parse()
            .map_err(|_| format!("invalid cluster number in UC record: {line}"))?;
        let record = |minus, cigar| {
            Ok(Some(ClusterRecord {
                cluster,
                label: fields[8],
                minus,
                cigar,
            }))
        };
        match fields[0] {
            "S" if cluster == clusters => record(false, ""),
            "H" if cluster < clusters => record(fields[4] == "-", fields[7]),
            "S" | "H" => Err(format!(
                "cluster {cluster} out of order in UC record: {line}"
            )),
            _ => Err(format!("not a cluster record: {line}")),
        }
    }
}
