use bio::io::fasta::Record;
use bio::io::fastq;
use std::collections::HashMap;
use std::sync::LazyLock;

use crate::fastx;
use crate::label::abundance;

#[derive(Clone, Debug, Default)]
/// A struct to store the size and annotation of sequences.
///
/// `merged` holds the label of every merged sequence, annotations
/// included, and its size is the sum of their `;size=` annotations. `seq`
/// and `qual` are those of the first merged sequence; `qual` is `None` for
/// FASTA input.
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

/// `error_probability` of every byte, computed once: merging looks up two
/// probabilities per position of every merged read.
static ERROR_PROBABILITY: LazyLock<[f64; 256]> =
    LazyLock::new(|| std::array::from_fn(|symbol| error_probability(symbol as u8)));

/// The quality symbol of an error probability: the Phred score
/// `trunc(-10 log10(p))`, limited to `MIN_QUALITY..=MAX_QUALITY`.
fn quality_symbol(probability: f64) -> u8 {
    let q = (-10.0 * probability.log10())
        .trunc()
        .clamp(f64::from(MIN_QUALITY), MAX_QUALITY);
    q as u8 + PHRED_OFFSET
}

/// Below this score a quality carries too little information to keep
/// against a better read.
const MIN_KEPT_QUALITY: u8 = 2;
/// A read that much less abundant than its cluster cannot raise a score.
const NEGLIGIBLE_SIZE_RATIO: u64 = 4;
/// Up to this cluster size, the abundance weights are exact in an `f64`.
const RESOLVABLE_SIZE: u64 = 1 << 40;

/// raises Q0 to Q1, as vsearch does on output.
fn floor_quality(qual: &[u8]) -> Vec<u8> {
    qual.iter()
        .map(|&q| q.max(PHRED_OFFSET + MIN_QUALITY))
        .collect()
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
        Cluster {
            size: abundance(record.id()),
            merged: vec![record.id().to_string()],
            seq: record.seq().to_vec(),
            qual: Some(floor_quality(record.qual())),
        }
    }

    /// merges `other` into this cluster: adds up the sizes, appends the
    /// merged names and combines the quality scores as the size-weighted
    /// average of the error probabilities. `reverse` tells that `other` is
    /// the reverse complement of this cluster's sequence.
    ///
    /// As in vsearch, a better score from a read at most a quarter of the
    /// cluster's size leaves an informative score unchanged. Such a read
    /// moves the average by less than one, and computing it could still
    /// lower the score by one through rounding.
    pub fn merge(&mut self, other: Cluster, reverse: bool) {
        self.merge_quality(other.size, other.qual.as_deref(), reverse);
        self.size += other.size;
        self.merged.extend(other.merged);
    }

    /// merges a read named `name` of abundance `size`, as `merge` does,
    /// without building a cluster for it first. `qual` is `None` for FASTA.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::cluster::Cluster;
    /// # use rsearch::fastx::Record;
    /// let mut cluster = Cluster::build_from_fastx(&Record::with_attrs("a", None, b"ACGT", b"I#5+"));
    /// cluster.merge_read("b", 1, Some(b"#I5,"), false);
    ///
    /// assert_eq!(cluster.size, 2);
    /// assert_eq!(cluster.merged, ["a", "b"]);
    /// assert_eq!(cluster.qual.as_deref(), Some(&b"&&5+"[..]));
    /// ```
    pub fn merge_read(&mut self, name: &str, size: u64, qual: Option<&[u8]>, reverse: bool) {
        self.merge_quality(size, qual, reverse);
        self.size += size;
        self.merged.push(name.to_string());
    }

    /// combines the qualities of a read of abundance `size` into the
    /// cluster's, before the sizes are added up.
    fn merge_quality(&mut self, size: u64, qual: Option<&[u8]>, reverse: bool) {
        if let (Some(q1), Some(q2)) = (&mut self.qual, qual) {
            let negligible =
                size <= self.size / NEGLIGIBLE_SIZE_RATIO && self.size < RESOLVABLE_SIZE;
            let (s1, s2) = (self.size as f64, size as f64);
            let probability = &*ERROR_PROBABILITY;
            let len = q1.len();
            for (i, q) in q1.iter_mut().enumerate() {
                let p = q2[if reverse { len - 1 - i } else { i }];
                if *q == p || (negligible && p > *q && *q >= PHRED_OFFSET + MIN_KEPT_QUALITY) {
                    continue;
                }
                let mixed = (probability[usize::from(*q)] * s1 + probability[usize::from(p)] * s2)
                    / (s1 + s2);
                *q = quality_symbol(mixed);
            }
        }
    }

    /// returns a Cluster from a `fastx::Record`, with qualities if the
    /// record is from FASTQ.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::cluster::Cluster;
    /// # use rsearch::fastx::Record;
    /// let cluster = Cluster::build_from_fastx(&Record::fasta_with_attrs("foo;size=3", None, b"ACGT"));
    ///
    /// assert_eq!(cluster.size, 3);
    /// assert_eq!(cluster.merged, ["foo;size=3"]);
    /// assert_eq!(cluster.qual, None);
    /// ```
    pub fn build_from_fastx(record: &fastx::Record) -> Cluster {
        Cluster {
            size: abundance(record.id()),
            merged: vec![record.id().to_string()],
            seq: record.seq().to_vec(),
            qual: record.is_fastq().then(|| floor_quality(record.qual())),
        }
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
    /// assert_eq!(cluster.merged[0], "foo;size=100");
    /// # }
    /// ```
    pub fn build_from_record(record: &Record) -> Cluster {
        Cluster {
            size: abundance(record.id()),
            merged: vec![record.id().to_string()],
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
    /// The normalized sequence being added, and its reverse complement.
    key: Vec<u8>,
    reverse: Vec<u8>,
    /// Leave out qualities, which are only needed to write FASTQ.
    no_qualities: bool,
}

impl Clusters {
    /// returns empty clusters that leave out the qualities of FASTQ
    /// records, skipping the work of merging them.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::cluster::Clusters;
    /// # use rsearch::fastx::Record;
    /// let mut clusters = Clusters::without_qualities();
    /// clusters.add_fastx(&Record::with_attrs("a", None, b"ACGT", b"IIII"), false);
    ///
    /// assert_eq!(clusters.into_clusters()[0].qual, None);
    /// ```
    pub fn without_qualities() -> Clusters {
        Clusters {
            no_qualities: true,
            ..Clusters::default()
        }
    }

    /// finds the cluster of `seq` by its normalized sequence, or with
    /// `strand`, its reverse complement; returns its index and whether it
    /// matched the reverse complement. Leaves the normalized sequence in
    /// `self.key`.
    fn find(&mut self, seq: &[u8], strand: bool) -> Option<(usize, bool)> {
        normalize_into(&mut self.key, seq);
        if let Some(&i) = self.index.get(self.key.as_slice()) {
            return Some((i, false));
        }
        if strand {
            reverse_complement_into(&mut self.reverse, &self.key);
            if let Some(&i) = self.index.get(self.reverse.as_slice()) {
                return Some((i, true));
            }
        }
        None
    }

    fn push(&mut self, cluster: Cluster) {
        self.index.insert(self.key.clone(), self.clusters.len());
        self.clusters.push(cluster);
    }

    /// adds a sequence, merging it into an existing cluster with the same
    /// normalized sequence, or with `strand`, its reverse complement.
    pub fn add(&mut self, cluster: Cluster, strand: bool) {
        match self.find(&cluster.seq, strand) {
            Some((i, reverse)) => self.clusters[i].merge(cluster, reverse),
            None => self.push(cluster),
        }
    }

    /// adds a record as `add` does, but builds a cluster only for a new
    /// sequence. A duplicate is merged without copying its sequence.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::cluster::Clusters;
    /// # use rsearch::fastx::Record;
    /// let mut clusters = Clusters::default();
    /// for (id, seq) in [("a", b"AACC"), ("b", b"aacc"), ("c", b"GGTT")] {
    ///     clusters.add_fastx(&Record::fasta_with_attrs(id, None, seq), true);
    /// }
    ///
    /// let clusters = clusters.into_clusters();
    /// assert_eq!(clusters.len(), 1);
    /// assert_eq!(clusters[0].merged, ["a", "b", "c"]);
    /// ```
    pub fn add_fastx(&mut self, record: &fastx::Record, strand: bool) {
        self.add_read(record.id(), abundance(record.id()), record, strand);
    }

    /// adds a record as `add_fastx` does, under the given label and
    /// abundance instead of its id and `;size=` annotation.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::cluster::Clusters;
    /// # use rsearch::fastx::Record;
    /// let mut clusters = Clusters::default();
    /// for header in ["a;size=3 one", "b;size=5 two"] {
    ///     let record = Record::fasta_with_attrs(header, None, b"ACGT");
    ///     clusters.add_read(record.header(), 1, &record, false);
    /// }
    ///
    /// let clusters = clusters.into_clusters();
    /// assert_eq!(clusters[0].size, 2);
    /// assert_eq!(clusters[0].merged, ["a;size=3 one", "b;size=5 two"]);
    /// ```
    pub fn add_read(&mut self, label: &str, size: u64, record: &fastx::Record, strand: bool) {
        let qual = (record.is_fastq() && !self.no_qualities).then(|| record.qual());
        match self.find(record.seq(), strand) {
            Some((i, reverse)) => self.clusters[i].merge_read(label, size, qual, reverse),
            None => self.push(Cluster {
                size,
                merged: vec![label.to_string()],
                seq: record.seq().to_vec(),
                qual: qual.map(floor_quality),
            }),
        }
    }

    /// returns the clusters in order of first appearance.
    pub fn into_clusters(self) -> Vec<Cluster> {
        self.clusters
    }
}

/// sorts clusters by decreasing size. The sort is stable, and ties are
/// broken by the label of the first merged sequence, then by first
/// appearance, as in vsearch.
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
    let mut normalized = Vec::with_capacity(seq.len());
    normalize_into(&mut normalized, seq);
    normalized
}

/// writes `normalize(seq)` into `dest`, reusing its allocation.
pub fn normalize_into(dest: &mut Vec<u8>, seq: &[u8]) {
    dest.clear();
    dest.extend(seq.iter().map(|&c| normalize_base(c)));
}

/// returns one normalized base: upper case, U as T, and N for anything
/// that is not an IUPAC nucleotide code.
pub fn normalize_base(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'U' => b'T',
        c @ (b'A' | b'B' | b'C' | b'D' | b'G' | b'H' | b'K' | b'M' | b'N' | b'R' | b'S' | b'T'
        | b'V' | b'W' | b'Y') => c,
        _ => b'N',
    }
}

/// returns the reverse complement of a normalized sequence.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::reverse_complement;
/// assert_eq!(reverse_complement(b"AACGTRN"), b"NYACGTT");
/// ```
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut complement = Vec::with_capacity(seq.len());
    reverse_complement_into(&mut complement, seq);
    complement
}

/// writes `reverse_complement(seq)` into `dest`, reusing its allocation.
pub fn reverse_complement_into(dest: &mut Vec<u8>, seq: &[u8]) {
    dest.clear();
    dest.extend(seq.iter().rev().map(|c| match c {
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
    }));
}

#[cfg(test)]
mod test {
    use super::*;

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
        assert_eq!(cluster.merged, ["foo;size=3"]);
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
        assert_eq!(a.merged, ["a;size=2", "b;size=5"]);
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
