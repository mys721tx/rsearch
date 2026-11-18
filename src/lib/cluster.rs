use bio::io::fasta::Record;
use bio::io::fastq;
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
/// A struct to store the size and annotation of sequences.
///
/// `seq` and `qual` are those of the first merged sequence; `qual` is
/// `None` for FASTA input.
pub struct Cluster {
    pub size: u64,
    pub merged: Vec<String>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

/// ASCII offset of Phred quality scores.
const PHRED_OFFSET: u8 = 33;
/// Lowest Phred score written: a score of 0 claims a certain error, so it
/// is raised to 1, as vsearch does.
const MIN_QUALITY: u8 = 1;
/// Highest Phred score representable in printable ASCII.
const MAX_QUALITY: f64 = 93.0;
/// Error probabilities are capped at a random guess among four bases.
const MAX_ERROR: f64 = 0.75;

fn error_probability(symbol: u8) -> f64 {
    10f64
        .powf(-f64::from(symbol.saturating_sub(PHRED_OFFSET)) / 10.0)
        .min(MAX_ERROR)
}

fn quality_symbol(probability: f64) -> u8 {
    let q = (-10.0 * probability.log10())
        .trunc()
        .clamp(f64::from(MIN_QUALITY), MAX_QUALITY);
    q as u8 + PHRED_OFFSET
}

impl Cluster {
    /// returns a Cluster from a `bio::io::fastq::Record`.
    ///
    /// # Example
    /// ```rust
    /// # use bio::io::fastq::Record;
    /// # use rsearch::cluster::Cluster;
    /// let record = Record::with_attrs("foo;size=3", None, b"ACGT", b"IIII");
    ///
    /// let cluster = Cluster::build_from_fastq_record(&record);
    ///
    /// assert_eq!(cluster.size, 3);
    /// assert_eq!(cluster.qual.as_deref(), Some(&b"IIII"[..]));
    /// ```
    pub fn build_from_fastq_record(record: &fastq::Record) -> Cluster {
        let (name, size) = parse_id(record.id());

        Cluster {
            size,
            merged: vec![name.to_string()],
            seq: record.seq().to_vec(),
            qual: Some(
                record
                    .qual()
                    .iter()
                    .map(|&q| q.max(PHRED_OFFSET + MIN_QUALITY))
                    .collect(),
            ),
        }
    }

    /// merges `other` into this cluster: adds up the sizes, appends the
    /// merged names and combines the quality scores as the size-weighted
    /// average of the error probabilities. `reverse` tells that `other` is
    /// the reverse complement of this cluster's sequence.
    pub fn merge(&mut self, other: Cluster, reverse: bool) {
        if let (Some(q1), Some(q2)) = (&mut self.qual, &other.qual) {
            let (s1, s2) = (self.size as f64, other.size as f64);
            let len = q1.len();
            for (i, q) in q1.iter_mut().enumerate() {
                let p = q2[if reverse { len - 1 - i } else { i }];
                if *q != p {
                    let mixed =
                        (error_probability(*q) * s1 + error_probability(p) * s2) / (s1 + s2);
                    *q = quality_symbol(mixed);
                }
            }
        }
        self.size += other.size;
        self.merged.extend(other.merged);
    }
    /// returns a Cluster from a `bio::io::fasta::Record`.
    ///
    /// # Example
    /// ```rust
    /// # extern crate bio;
    /// # extern crate rsearch;
    /// # use bio::io::fasta::Record;
    /// # use rsearch::cluster::Cluster;
    /// # fn main() {
    /// let record = Record::with_attrs("foo;size=100", Option::None, b"AAAA");
    ///
    /// let cluster = Cluster::build_from_record(&record);
    ///
    /// assert_eq!(cluster.size, 100);
    /// assert_eq!(cluster.merged[0], "foo");
    /// # }
    /// ```
    pub fn build_from_record(record: &Record) -> Cluster {
        let (name, size) = parse_id(record.id());

        let merged = vec![name.to_string()];
        Cluster {
            size,
            merged,
            seq: record.seq().to_vec(),
            qual: None,
        }
    }
}

/// Dereplicated sequences: the clusters in order of first appearance, and
/// the normalized sequence of each cluster mapped to its index.
///
/// # Example
/// ```rust
/// # use bio::io::fasta::Record;
/// # use rsearch::cluster::{Cluster, Clusters};
/// let mut clusters = Clusters::default();
/// for (id, seq) in [("a", b"AACC"), ("b", b"aacc"), ("c", b"GGTT")] {
///     clusters.add(Cluster::build_from_record(&Record::with_attrs(id, None, seq)), true);
/// }
///
/// let clusters = clusters.into_clusters();
/// assert_eq!(clusters.len(), 1);
/// assert_eq!(clusters[0].size, 3);
/// ```
#[derive(Debug, Default)]
pub struct Clusters {
    clusters: Vec<Cluster>,
    index: HashMap<Vec<u8>, usize>,
}

impl Clusters {
    /// adds a sequence, merging it into an existing cluster with the same
    /// normalized sequence, or with `strand`, its reverse complement.
    pub fn add(&mut self, cluster: Cluster, strand: bool) {
        let key = normalize(&cluster.seq);

        if let Some(&i) = self.index.get(&key) {
            self.clusters[i].merge(cluster, false);
            return;
        }

        if strand && let Some(&i) = self.index.get(&reverse_complement(&key)) {
            self.clusters[i].merge(cluster, true);
            return;
        }

        self.index.insert(key, self.clusters.len());
        self.clusters.push(cluster);
    }

    /// returns the clusters in order of first appearance.
    pub fn into_clusters(self) -> Vec<Cluster> {
        self.clusters
    }
}

/// sorts clusters by decreasing size. The sort is stable, and ties are
/// broken by name, then by first appearance.
pub fn sort_by_size(clusters: &mut [Cluster]) {
    clusters.sort_by(|a, b| {
        b.size
            .cmp(&a.size)
            .then_with(|| a.merged[0].cmp(&b.merged[0]))
    });
}

/// returns the sequence in upper case with U replaced by T and every
/// symbol that is not an IUPAC nucleotide code replaced by N, so that
/// sequences differing only in case or U/T compare equal.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::normalize;
/// assert_eq!(normalize(b"acgUrx"), b"ACGTRN");
/// ```
pub fn normalize(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|c| match c.to_ascii_uppercase() {
            b'U' => b'T',
            c @ (b'A' | b'B' | b'C' | b'D' | b'G' | b'H' | b'K' | b'M' | b'N' | b'R' | b'S'
            | b'T' | b'V' | b'W' | b'Y') => c,
            _ => b'N',
        })
        .collect()
}

/// returns the reverse complement of a normalized sequence.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::reverse_complement;
/// assert_eq!(reverse_complement(b"AACGTRN"), b"NYACGTT");
/// ```
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'S' => b'S',
            b'W' => b'W',
            _ => b'N',
        })
        .collect()
}

/// returns the name and the size of a semicolon-saperated string.
///
/// `parse_id` parses a semicolon-saperated string annotated sequence ID and
/// returns the name and the size of that sequence.
///
/// `parse_id` takes the first non-empty item without an equal sign as the name
/// of the sequence. `parse_id` treat non-empty item with a single equal sign as
/// a key-value pair. `parse_id` takes the value of the first pair with "size" as
/// key as the size of the sequence.
///
/// If no valid item exists, the name defaults to "sequence" and the size
/// defaults to 1.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::parse_id;
/// # fn main() {
/// let (name, size) = parse_id("foo;size=100");
/// assert_eq!(name, "foo");
/// assert_eq!(size, 100);
/// # }
/// ```
pub fn parse_id(id: &str) -> (&str, u64) {
    let mut name = None;
    let mut size: Option<u64> = None;

    for item in id.split(';') {
        let pair: Vec<&str> = item.split('=').collect();
        match pair.len() {
            1 => {
                if name.is_none() && !pair[0].is_empty() {
                    name = Some(pair[0]);
                }
            }
            2 => {
                if size.is_none()
                    && pair[0] == "size"
                    && let Ok(x) = pair[1].parse().map(|x| if x > 0 { x } else { 1 })
                {
                    size = Some(x);
                }
            }
            _ => (),
        }
    }

    (name.unwrap_or("sequence"), size.unwrap_or(1))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_id() {
        let mut input = "foo;size=100";
        let mut expected = ("foo", 100);
        assert_eq!(parse_id(input), expected,);

        input = "species=9606;spam=fried;size=23";
        expected = ("sequence", 23);
        assert_eq!(parse_id(input), expected,);

        input = "egg=scrumbled;spam=fried;pancake;size=127;waffle";
        expected = ("pancake", 127);
        assert_eq!(parse_id(input), expected,);

        input = "かえで;lily;Ἶρῐς";
        expected = ("かえで", 1);
        assert_eq!(parse_id(input), expected,);

        input = "size=0022;size=2;size=3;size=4";
        expected = ("sequence", 22);
        assert_eq!(parse_id(input), expected,);

        input = "size=abc;size=甲乙丙;size=αβγ";
        expected = ("sequence", 1);
        assert_eq!(parse_id(input), expected,);

        input = "size;size=200;size=400=300;";
        expected = ("size", 200);
        assert_eq!(parse_id(input), expected,);
    }

    fn fastq_cluster(id: &str, seq: &[u8], qual: &[u8]) -> Cluster {
        Cluster::build_from_fastq_record(&fastq::Record::with_attrs(id, None, seq, qual))
    }

    #[test]
    fn test_normalize() {
        // case and U/T are ignored
        assert_eq!(normalize(b"acgtu"), b"ACGTT");
        assert_eq!(normalize(b"ACGU"), normalize(b"acgt"));
        // all IUPAC codes are kept, anything else becomes N
        assert_eq!(normalize(b"RYSWKMBDHVN"), b"RYSWKMBDHVN");
        assert_eq!(normalize(b"rx*-"), b"RNNN");
        // ambiguity codes stay distinct from each other and from N
        assert_ne!(normalize(b"ACGTR"), normalize(b"ACGTN"));
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"AACC"), b"GGTT");
        assert_eq!(reverse_complement(b"ACGTRYKMBVDHSWN"), b"NWSDHBVKMRYACGT");
        // complementing twice gives back the sequence
        let seq = normalize(b"ACGTRYKMBVDHSWN");
        assert_eq!(reverse_complement(&reverse_complement(&seq)), seq);
    }

    #[test]
    fn test_build_from_record() {
        let record = Record::with_attrs("foo;size=3", Some("desc"), b"acgu");
        let cluster = Cluster::build_from_record(&record);
        assert_eq!(cluster.size, 3);
        assert_eq!(cluster.merged, ["foo"]);
        // the sequence is kept as written; matching uses `normalize`
        assert_eq!(cluster.seq, b"acgu");
        assert_eq!(cluster.qual, None);
    }

    #[test]
    fn test_build_from_fastq_record() {
        let cluster = fastq_cluster("bar", b"ACGT", b"!\"+I");
        assert_eq!(cluster.size, 1);
        assert_eq!(cluster.merged, ["bar"]);
        // Q0 is raised to Q1, as vsearch does
        assert_eq!(cluster.qual.as_deref(), Some(&b"\"\"+I"[..]));
    }

    #[test]
    fn test_merge_fasta() {
        let mut a = Cluster::build_from_record(&Record::with_attrs("a;size=2", None, b"ACGT"));
        let b = Cluster::build_from_record(&Record::with_attrs("b;size=5", None, b"acgt"));
        a.merge(b, false);
        assert_eq!(a.size, 7);
        assert_eq!(a.merged, ["a", "b"]);
        // the first sequence represents the cluster
        assert_eq!(a.seq, b"ACGT");
        assert_eq!(a.qual, None);
    }

    // Expected qualities below are the output of vsearch --fastx_uniques.

    #[test]
    fn test_merge_quality_equal_sizes() {
        let mut a = fastq_cluster("a", b"ACGT", b"I#5+");
        a.merge(fastq_cluster("b", b"ACGT", b"#I5,"), false);
        assert_eq!(a.qual.as_deref(), Some(&b"&&5+"[..]));
    }

    #[test]
    fn test_merge_quality_weighted_by_size() {
        let mut a = fastq_cluster("a;size=3", b"ACGT", b"5555");
        a.merge(fastq_cluster("b;size=2", b"acgu", b"++++"), false);
        a.merge(fastq_cluster("c", b"ACGT", b"????"), false);
        assert_eq!(a.size, 6);
        assert_eq!(a.qual.as_deref(), Some(&b"...."[..]));

        // a single better read does not lift a much more abundant cluster
        let mut a = fastq_cluster("a;size=4", b"ACGT", b"++++");
        a.merge(fastq_cluster("b", b"ACGT", b"IIII"), false);
        assert_eq!(a.qual.as_deref(), Some(&b"++++"[..]));
    }

    #[test]
    fn test_merge_quality_reverse_strand() {
        // the qualities of a reverse-complemented read are reversed
        let mut a = fastq_cluster("a", b"AACC", b"ABCD");
        a.merge(fastq_cluster("b", b"GGTT", b"IIII"), true);
        assert_eq!(a.qual.as_deref(), Some(&b"CDEE"[..]));

        let mut a = fastq_cluster("a", b"AACC", b"IIII");
        a.merge(fastq_cluster("b", b"GGTT", b"#+5?"), true);
        let mut b = fastq_cluster("a", b"AACC", b"IIII");
        b.merge(fastq_cluster("b", b"AACC", b"?5+#"), false);
        assert_eq!(a.qual, b.qual);
    }

    #[test]
    fn test_merge_quality_bounds() {
        // error probabilities are capped at 0.75 and scores floored at Q1
        let mut a = fastq_cluster("a", b"AC", b"!#");
        a.merge(fastq_cluster("b", b"AC", b"#!"), false);
        assert_eq!(a.qual.as_deref(), Some(&b"\"\""[..]));

        // the highest score is kept when identical
        let mut a = fastq_cluster("a", b"A", b"~");
        a.merge(fastq_cluster("b", b"A", b"~"), false);
        assert_eq!(a.qual.as_deref(), Some(&b"~"[..]));
    }
}
