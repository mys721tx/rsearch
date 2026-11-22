use bio::io::fasta::Record;
use bio::io::fastq;
use std::collections::{BTreeMap, HashMap};
use std::hash::BuildHasherDefault;
use std::sync::LazyLock;

use crate::fastx;
use crate::label::abundance;

#[derive(Clone, Debug, Default)]
/// A struct to store the size and annotation of sequences.
///
/// `merged` holds the label of every merged sequence, annotations
/// included, and its size is the sum of their abundances. `reverse` and
/// `lengths` run parallel to it: whether each sequence was merged through
/// its reverse complement, and its length. `seq` and `qual` are those of
/// the first merged sequence; `qual` is `None` for FASTA input, and keeps
/// the input encoding until `Quality::encode` writes it out.
pub struct Cluster {
    pub size: u64,
    pub merged: Vec<String>,
    pub reverse: Vec<bool>,
    pub lengths: Vec<usize>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

/// Error probabilities are capped at a random guess among four bases.
const MAX_ERROR: f64 = 0.75;
/// Lowest score kept. A score of 0 claims a certain error and is raised to
/// 1, as vsearch does.
const MIN_QUALITY: u8 = 1;
/// Below this score a quality carries too little information to keep
/// against a better read.
const MIN_KEPT_QUALITY: i16 = 2;
/// A read that much less abundant than its cluster cannot raise a score.
const NEGLIGIBLE_SIZE_RATIO: u64 = 4;
/// Up to this cluster size, the abundance weights are exact in an `f64`.
const RESOLVABLE_SIZE: u64 = 1 << 40;

/// How quality scores are read, merged and written.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::Quality;
/// // Phred+64 in, Phred+33 out, scores limited to 2..=40
/// let quality = Quality {
///     offset: 64,
///     offset_out: 33,
///     min_out: 2,
///     max_out: 40,
///     ..Quality::default()
/// };
/// assert_eq!(quality.encode(b"@Ahi"), b"##II");
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Quality {
    /// ASCII offset of the input scores.
    pub offset: u8,
    /// ASCII offset of the output scores.
    pub offset_out: u8,
    /// Lowest score written.
    pub min_out: u8,
    /// Highest score written.
    pub max_out: u8,
    /// Keep the best score of the merged sequences instead of the
    /// size-weighted average of their error probabilities.
    pub keep_best: bool,
}

impl Default for Quality {
    /// Phred+33, with every score the encoding can carry.
    fn default() -> Quality {
        Quality {
            offset: 33,
            offset_out: 33,
            min_out: 0,
            max_out: 93,
            keep_best: false,
        }
    }
}

impl Quality {
    /// the score of a symbol, which may be negative.
    fn score(&self, symbol: u8) -> i16 {
        i16::from(symbol) - i16::from(self.offset)
    }

    /// the error probability of every symbol, capped at `MAX_ERROR`.
    fn probabilities(&self) -> [f64; 256] {
        std::array::from_fn(|symbol| {
            10f64
                .powf(-f64::from(self.score(symbol as u8)) / 10.0)
                .min(MAX_ERROR)
        })
    }

    /// the input symbol of an error probability: the score
    /// `trunc(-10 log10(p))`, limited to `min_out..=max_out`.
    fn symbol(&self, probability: f64) -> u8 {
        let q = (-10.0 * probability.log10()).trunc();
        let q = q.min(f64::from(self.max_out)).max(f64::from(self.min_out));
        (q as u8).saturating_add(self.offset)
    }

    /// the output symbol of every stored symbol: its score raised to at
    /// least 1, limited to `min_out..=max_out` and written with
    /// `offset_out`.
    pub fn output_table(&self) -> [u8; 256] {
        std::array::from_fn(|symbol| {
            let q = self
                .score(symbol as u8)
                .max(i16::from(MIN_QUALITY))
                .min(i16::from(self.max_out))
                .max(i16::from(self.min_out));
            (q + i16::from(self.offset_out)).clamp(0, 255) as u8
        })
    }

    /// re-encodes stored scores for output, as `output_table` does.
    pub fn encode(&self, qual: &[u8]) -> Vec<u8> {
        let table = self.output_table();
        qual.iter().map(|&q| table[usize::from(q)]).collect()
    }
}

/// The error probabilities of the default encoding.
static DEFAULT_PROBABILITY: LazyLock<[f64; 256]> =
    LazyLock::new(|| Quality::default().probabilities());

impl Cluster {
    /// returns a cluster of one sequence, without qualities.
    fn new(label: &str, size: u64, seq: &[u8]) -> Cluster {
        Cluster {
            size,
            merged: vec![label.to_string()],
            reverse: vec![false],
            lengths: vec![seq.len()],
            seq: seq.to_vec(),
            qual: None,
        }
    }

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
            qual: Some(record.qual().to_vec()),
            ..Cluster::new(record.id(), abundance(record.id()), record.seq())
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
        let quality = Quality::default();
        self.merge_quality(
            other.size,
            other.qual.as_deref(),
            reverse,
            &quality,
            &DEFAULT_PROBABILITY,
        );
        self.size += other.size;
        self.merged.extend(other.merged);
        self.reverse
            .extend(other.reverse.iter().map(|&r| r != reverse));
        self.lengths.extend(other.lengths);
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
        let quality = Quality::default();
        self.merge_with(name, size, qual, reverse, &quality, &DEFAULT_PROBABILITY);
    }

    /// `merge_read` with the given quality settings and their error
    /// probabilities.
    fn merge_with(
        &mut self,
        name: &str,
        size: u64,
        qual: Option<&[u8]>,
        reverse: bool,
        quality: &Quality,
        probability: &[f64; 256],
    ) {
        self.merge_quality(size, qual, reverse, quality, probability);
        self.size += size;
        self.merged.push(name.to_string());
        self.reverse.push(reverse);
        self.lengths.push(self.seq.len());
    }

    /// combines the qualities of a read of abundance `size` into the
    /// cluster's, before the sizes are added up.
    fn merge_quality(
        &mut self,
        size: u64,
        qual: Option<&[u8]>,
        reverse: bool,
        quality: &Quality,
        probability: &[f64; 256],
    ) {
        if let (Some(q1), Some(q2)) = (&mut self.qual, qual) {
            let len = q1.len();
            let other = |i: usize| q2[if reverse { len - 1 - i } else { i }];
            if quality.keep_best {
                for (i, q) in q1.iter_mut().enumerate() {
                    *q = (*q).max(other(i));
                }
                return;
            }
            let negligible =
                size <= self.size / NEGLIGIBLE_SIZE_RATIO && self.size < RESOLVABLE_SIZE;
            // the stored symbols a negligible read cannot raise: informative
            // scores that need no clamping on output
            let lowest = (i16::from(quality.offset)
                + MIN_KEPT_QUALITY.max(i16::from(quality.min_out)))
            .clamp(0, 256) as u16;
            let highest = u16::from(quality.offset) + u16::from(quality.max_out);
            let (s1, s2) = (self.size as f64, size as f64);
            for (i, q) in q1.iter_mut().enumerate() {
                let p = other(i);
                if *q == p || (negligible && p > *q && (lowest..=highest).contains(&u16::from(*q)))
                {
                    continue;
                }
                let mixed = (probability[usize::from(*q)] * s1 + probability[usize::from(p)] * s2)
                    / (s1 + s2);
                *q = quality.symbol(mixed);
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
            qual: record.is_fastq().then(|| record.qual().to_vec()),
            ..Cluster::new(record.id(), abundance(record.id()), record.seq())
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
        Cluster::new(record.id(), abundance(record.id()), record.seq())
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
#[derive(Debug)]
pub struct Clusters {
    clusters: Vec<Cluster>,
    index: HashMap<Vec<u8>, usize>,
    /// The key of the sequence being added (its normalized sequence, and
    /// with `by_label` its label), and that of its reverse complement.
    key: Vec<u8>,
    reverse: Vec<u8>,
    /// Leave out qualities, which are only needed to write FASTQ.
    no_qualities: bool,
    /// Merge only sequences with identical labels.
    by_label: bool,
    quality: Quality,
    probability: Box<[f64; 256]>,
}

impl Default for Clusters {
    fn default() -> Clusters {
        Clusters {
            clusters: Vec::new(),
            index: HashMap::new(),
            key: Vec::new(),
            reverse: Vec::new(),
            no_qualities: false,
            by_label: false,
            quality: Quality::default(),
            probability: Box::new(*DEFAULT_PROBABILITY),
        }
    }
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

    /// reads, merges and keeps qualities with these settings.
    pub fn with_quality(self, quality: Quality) -> Clusters {
        Clusters {
            quality,
            probability: Box::new(quality.probabilities()),
            ..self
        }
    }

    /// merges only sequences that also have identical labels, as vsearch
    /// --derep_id does.
    ///
    /// # Example
    /// ```rust
    /// # use rsearch::cluster::Clusters;
    /// # use rsearch::fastx::Record;
    /// let mut clusters = Clusters::default().by_label();
    /// for id in ["a", "b", "a"] {
    ///     clusters.add_fastx(&Record::fasta_with_attrs(id, None, b"ACGT"), false);
    /// }
    ///
    /// let sizes: Vec<u64> = clusters.into_clusters().iter().map(|c| c.size).collect();
    /// assert_eq!(sizes, [2, 1]);
    /// ```
    pub fn by_label(self) -> Clusters {
        Clusters {
            by_label: true,
            ..self
        }
    }

    /// finds the cluster of `seq` by its normalized sequence, or with
    /// `strand`, its reverse complement; with `by_label`, the label must
    /// match too. Returns its index and whether it matched the reverse
    /// complement. Leaves the key of the sequence in `self.key`.
    fn find(&mut self, seq: &[u8], label: &str, strand: bool) -> Option<(usize, bool)> {
        normalize_into(&mut self.key, seq);
        let sequence_len = self.key.len();
        if self.by_label {
            // a normalized sequence has no 0 byte and ends at the first
            self.key.push(0);
            self.key.extend_from_slice(label.as_bytes());
        }
        if let Some(&i) = self.index.get(self.key.as_slice()) {
            return Some((i, false));
        }
        if strand {
            reverse_complement_into(&mut self.reverse, &self.key[..sequence_len]);
            self.reverse.extend_from_slice(&self.key[sequence_len..]);
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
        match self.find(&cluster.seq, &cluster.merged[0], strand) {
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
        match self.find(record.seq(), label, strand) {
            Some((i, reverse)) => self.clusters[i].merge_with(
                label,
                size,
                qual,
                reverse,
                &self.quality,
                &self.probability,
            ),
            None => {
                let cluster = Cluster {
                    qual: qual.map(<[u8]>::to_vec),
                    ..Cluster::new(label, size, record.seq())
                };
                self.push(cluster);
            }
        }
    }

    /// returns the clusters in order of first appearance.
    pub fn into_clusters(self) -> Vec<Cluster> {
        self.clusters
    }
}

/// A sequence to dereplicate by prefix: its label, its abundance, the
/// abundance annotated in its label (which orders sequences of the same
/// length) and its sequence.
#[derive(Clone, Debug)]
pub struct PrefixRead {
    pub label: String,
    pub size: u64,
    pub annotated: u64,
    pub seq: Vec<u8>,
}

/// dereplicates sequences by prefix, as vsearch --derep_prefix does.
///
/// Sequences are taken from the shortest, then by decreasing annotated
/// abundance, label and input order. Each is merged into the cluster whose
/// representative is identical to it or, failing that, is its longest
/// prefix. In the latter case the longer sequence becomes the
/// representative and comes first among the members. Returns the clusters
/// in the order their representatives were taken.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::{PrefixRead, dereplicate_prefixes};
/// let read = |label: &str, seq: &[u8]| PrefixRead {
///     label: label.into(),
///     size: 1,
///     annotated: 1,
///     seq: seq.to_vec(),
/// };
/// let clusters = dereplicate_prefixes(vec![
///     read("a", b"ACGTAC"),
///     read("b", b"ACG"),
///     read("c", b"acgt"),
///     read("d", b"TTT"),
/// ]);
///
/// assert_eq!(clusters.len(), 2);
/// assert_eq!(clusters[1].seq, b"ACGTAC");
/// assert_eq!(clusters[1].merged, ["a", "c", "b"]);
/// assert_eq!(clusters[1].lengths, [6, 4, 3]);
/// ```
pub fn dereplicate_prefixes(reads: Vec<PrefixRead>) -> Vec<Cluster> {
    let mut order: Vec<usize> = (0..reads.len()).collect();
    order.sort_by(|&a, &b| {
        let (a, b) = (&reads[a], &reads[b]);
        a.seq
            .len()
            .cmp(&b.seq.len())
            .then(b.annotated.cmp(&a.annotated))
            .then_with(|| a.label.cmp(&b.label))
    });

    let mut reads: Vec<Option<PrefixRead>> = reads.into_iter().map(Some).collect();
    let mut clusters: Vec<Option<Cluster>> = Vec::new();
    // the normalized representative of each cluster; the live clusters by
    // the hash of their representative, as vsearch finds them; and how
    // many live representatives have each length
    let mut keys: Vec<Vec<u8>> = Vec::new();
    let mut index: HashMap<u64, Vec<usize>, BuildHasherDefault<PassThrough>> =
        HashMap::default();
    let mut lengths: BTreeMap<usize, usize> = BTreeMap::new();
    // the hash of each prefix of the read, the empty one first
    let mut hashes: Vec<u64> = Vec::new();

    for i in order {
        let read = reads[i].take().unwrap();
        let key = normalize(&read.seq);
        let len = key.len();
        prefix_hashes(&mut hashes, &key);
        // the live cluster whose representative is `key[..l]`
        let find = |l: usize| {
            index
                .get(&hashes[l])
                .and_then(|list| list.iter().copied().find(|&j| keys[j] == key[..l]))
        };

        if let Some(j) = find(len) {
            let cluster = clusters[j].as_mut().unwrap();
            cluster.size += read.size;
            cluster.merged.push(read.label);
            cluster.reverse.push(false);
            cluster.lengths.push(len);
            continue;
        }

        let prefix = lengths
            .range(..len)
            .rev()
            .find_map(|(&l, _)| find(l).map(|j| (l, j)));

        let mut cluster = Cluster::new(&read.label, read.size, &read.seq);
        if let Some((l, j)) = prefix {
            let list = index.get_mut(&hashes[l]).unwrap();
            list.retain(|&k| k != j);
            if list.is_empty() {
                index.remove(&hashes[l]);
            }
            if let Some(count) = lengths.get_mut(&l) {
                *count -= 1;
                if *count == 0 {
                    lengths.remove(&l);
                }
            }
            let old = clusters[j].take().unwrap();
            cluster.size += old.size;
            cluster.merged.extend(old.merged);
            cluster.reverse.extend(old.reverse);
            cluster.lengths.extend(old.lengths);
        }
        index.entry(hashes[len]).or_default().push(clusters.len());
        keys.push(key);
        *lengths.entry(len).or_insert(0) += 1;
        clusters.push(Some(cluster));
    }

    clusters.into_iter().flatten().collect()
}

/// writes the 64-bit FNV-1a hash of each prefix of `seq` to `hashes`,
/// from the empty prefix to the whole, in one pass, as vsearch's
/// `--derep_prefix` does.
fn prefix_hashes(hashes: &mut Vec<u64>, seq: &[u8]) {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    hashes.clear();
    let mut hash = OFFSET;
    hashes.push(hash);
    for &c in seq {
        hash = (hash ^ u64::from(c)).wrapping_mul(PRIME);
        hashes.push(hash);
    }
}

/// A hasher for keys that are hashes already. It keeps them as they are.
#[derive(Default)]
struct PassThrough(u64);

impl std::hash::Hasher for PassThrough {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 << 8) | u64::from(b);
        }
    }

    fn write_u64(&mut self, x: u64) {
        self.0 = x;
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

/// The complement of every byte, keeping its case, as vsearch's. U is
/// complemented to A, the IUPAC codes swap in pairs, and S, W and N are
/// their own complements. Bytes the reader does not let through map to N.
const COMPLEMENT: [u8; 256] = {
    let mut table = [b'N'; 256];
    let pairs: [(u8, u8); 10] = [
        (b'A', b'T'),
        (b'C', b'G'),
        (b'R', b'Y'),
        (b'K', b'M'),
        (b'B', b'V'),
        (b'D', b'H'),
        (b'S', b'S'),
        (b'W', b'W'),
        (b'N', b'N'),
        (b'U', b'A'),
    ];
    let mut i = 0;
    while i < pairs.len() {
        let (a, b) = pairs[i];
        table[a as usize] = b;
        table[a.to_ascii_lowercase() as usize] = b.to_ascii_lowercase();
        // the pairs swap both ways, but U only ever becomes A
        if a != b'U' {
            table[b as usize] = a;
            table[b.to_ascii_lowercase() as usize] = a.to_ascii_lowercase();
        }
        i += 1;
    }
    table
};

/// writes the reverse complement of `seq` into `dest`, keeping the case of
/// each base, as vsearch's `--fastx_revcomp` does.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::reverse_complement_keeping_case;
/// let mut dest = Vec::new();
/// reverse_complement_keeping_case(&mut dest, b"AcgURyn");
/// assert_eq!(dest, b"nrYAcgT");
/// ```
pub fn reverse_complement_keeping_case(dest: &mut Vec<u8>, seq: &[u8]) {
    dest.clear();
    dest.extend(seq.iter().rev().map(|&c| COMPLEMENT[usize::from(c)]));
}

/// The complement of a base, keeping its case, as vsearch's.
///
/// # Example
/// ```rust
/// # use rsearch::cluster::complement;
/// assert_eq!(complement(b'U'), b'A');
/// assert_eq!(complement(b'r'), b'y');
/// assert_eq!(complement(b'X'), b'N');
/// ```
pub fn complement(base: u8) -> u8 {
    COMPLEMENT[usize::from(base)]
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_complement_keeping_case() {
        let mut dest = Vec::new();
        reverse_complement_keeping_case(&mut dest, b"ACGUuRYKMBVDHSWNacgtrykmbvdhswn");
        assert_eq!(dest, b"nwsdhbvkmryacgtNWSDHBVKMRYaACGT");
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
        // qualities are kept as read, and Q0 is raised to Q1 on output, as
        // vsearch does
        assert_eq!(cluster.qual.as_deref(), Some(&b"!\"+I"[..]));
        let encoded = Quality::default().encode(cluster.qual.as_deref().unwrap());
        assert_eq!(encoded, b"\"\"+I");
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
        // error probabilities are capped at 0.75, and Q0 merges as Q1
        let mut a = fastq_cluster("a", b"AC", b"!#");
        a.merge(fastq_cluster("b", b"AC", b"#!"), false);
        assert_eq!(a.qual.as_deref(), Some(&b"\"\""[..]));

        // the highest score is kept when identical
        let mut a = fastq_cluster("a", b"A", b"~");
        a.merge(fastq_cluster("b", b"A", b"~"), false);
        assert_eq!(a.qual.as_deref(), Some(&b"~"[..]));
    }
}
