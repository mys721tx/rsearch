//! vsearch's multiple alignment of a cluster (`msa.cpp`). Each member's
//! alignment to the centroid is gapped by the longest insertion of any
//! member before each centroid position. The profile of each column is
//! weighted by abundance and gives the consensus.

use crate::alignment::parse_cigar;

/// The counters of a profile column: A, C, G, T or U, other IUPAC codes,
/// and gaps.
const A: usize = 0;
const T: usize = 3;
const N: usize = 4;
const GAP: usize = 5;

/// A member of a cluster: its sequence as aligned (reverse-complemented
/// for the minus strand), its weight in the profile, and its alignment to
/// the centroid as (operation, run) pairs: M both, D the member alone, I
/// the centroid alone.
pub struct Member {
    pub seq: Vec<u8>,
    pub weight: u64,
    pub ops: Vec<(u8, usize)>,
}

/// parses the CIGAR of a member of `member_len` aligned to a centroid of
/// `centroid_len`: `=` is a perfect alignment, and runs of one operation
/// are merged. Fails when it does not fit the two sequences.
///
/// # Example
/// ```rust
/// # use rsearch::msa::member_ops;
/// assert_eq!(member_ops("=", 4, 4).unwrap(), vec![(b'M', 4)]);
/// assert_eq!(member_ops("2MD2D2M", 7, 4).unwrap(), vec![(b'M', 2), (b'D', 3), (b'M', 2)]);
/// assert!(member_ops("=", 5, 4).is_err());
/// assert!(member_ops("4M", 4, 5).is_err());
/// ```
pub fn member_ops(
    cigar: &str,
    member_len: usize,
    centroid_len: usize,
) -> Result<Vec<(u8, usize)>, String> {
    if cigar == "=" {
        if member_len != centroid_len {
            return Err("a perfect alignment between sequences of different lengths".to_string());
        }
        return Ok(if member_len == 0 {
            Vec::new()
        } else {
            vec![(b'M', member_len)]
        });
    }
    let mut ops: Vec<(u8, usize)> = Vec::new();
    for (op, run) in parse_cigar(cigar)? {
        match ops.last_mut() {
            Some((last, total)) if *last == op => *total += run,
            _ => ops.push((op, run)),
        }
    }
    let (mut centroid, mut member) = (0usize, 0usize);
    for &(op, run) in &ops {
        if op != b'D' {
            centroid += run;
        }
        if op != b'I' {
            member += run;
        }
    }
    if (centroid, member) != (centroid_len, member_len) {
        return Err(format!(
            "CIGAR {cigar} does not fit sequences of {member_len} and {centroid_len}"
        ));
    }
    Ok(ops)
}

/// A cluster's alignment: a row a sequence, centroid first, and the
/// consensus row (`+` over the terminal insertions, `-` where gaps win),
/// the consensus sequence and the profile.
pub struct Alignment {
    pub rows: Vec<Vec<u8>>,
    pub consensus_row: Vec<u8>,
    pub consensus: Vec<u8>,
    /// for each column: A, C, G, T or U, other IUPAC codes, gaps.
    pub profile: Vec<[u64; 6]>,
}

impl Alignment {
    /// the profile line of column `i`: its position, consensus symbol and
    /// counts of A, C, G, T, gaps and other codes, as vsearch writes them.
    pub fn profile_line(&self, i: usize) -> String {
        let c = &self.profile[i];
        format!(
            "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.consensus_row[i] as char, c[0], c[1], c[2], c[3], c[GAP], c[N]
        )
    }
}

/// A row being written, and the profile it adds to.
struct Row<'a> {
    row: Vec<u8>,
    profile: &'a mut [[u64; 6]],
    weight: u64,
}

impl Row<'_> {
    fn put(&mut self, c: u8) {
        let column = self.row.len();
        let counter = match c.to_ascii_uppercase() {
            b'A' => Some(A),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(T),
            b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N' => Some(N),
            b'-' => Some(GAP),
            _ => None,
        };
        if let Some(k) = counter {
            self.profile[column][k] += self.weight;
        }
        self.row.push(c);
    }

    fn gaps(&mut self, n: usize) {
        for _ in 0..n {
            self.put(b'-');
        }
    }
}

/// aligns a cluster: the centroid of `weight`, and its members.
pub fn align(centroid: &[u8], weight: u64, members: &[Member]) -> Alignment {
    let len = centroid.len();
    // the longest insertion before each centroid position, and after it
    let mut insertions = vec![0usize; len + 1];
    for m in members {
        let mut pos = 0;
        for &(op, run) in &m.ops {
            if op == b'D' {
                insertions[pos] = insertions[pos].max(run);
            } else {
                pos += run;
            }
        }
    }
    let width = len + insertions.iter().sum::<usize>();
    let mut profile = vec![[0u64; 6]; width];
    let mut rows = Vec::with_capacity(members.len() + 1);

    let mut row = Row {
        row: Vec::with_capacity(width),
        profile: &mut profile,
        weight,
    };
    for (i, &c) in centroid.iter().enumerate() {
        row.gaps(insertions[i]);
        row.put(c);
    }
    row.gaps(insertions[len]);
    rows.push(row.row);

    for m in members {
        let mut row = Row {
            row: Vec::with_capacity(width),
            profile: &mut profile,
            weight: m.weight,
        };
        let (mut qpos, mut tpos) = (0usize, 0usize);
        let mut inserted = false;
        for &(op, run) in &m.ops {
            match op {
                b'D' => {
                    for _ in 0..run {
                        row.put(m.seq[tpos]);
                        tpos += 1;
                    }
                    row.gaps(insertions[qpos].saturating_sub(run));
                    inserted = true;
                }
                b'M' => {
                    for _ in 0..run {
                        if !inserted {
                            row.gaps(insertions[qpos]);
                        }
                        row.put(m.seq[tpos]);
                        tpos += 1;
                        qpos += 1;
                        inserted = false;
                    }
                }
                _ => {
                    for _ in 0..run {
                        if !inserted {
                            row.gaps(insertions[qpos]);
                        }
                        row.put(b'-');
                        qpos += 1;
                        inserted = false;
                    }
                }
            }
        }
        if !inserted {
            row.gaps(insertions[qpos]);
        }
        rows.push(row.row);
    }

    // the consensus of the columns between the terminal insertions
    const SYMBOLS: &[u8; 16] = b"-ACMGRSVTWYHKDBN";
    let (left, right) = (insertions[0], insertions[len]);
    let mut consensus_row = vec![b'+'; width];
    let mut consensus = Vec::with_capacity(width);
    for i in left..width - right {
        let counts = &profile[i];
        // the first of the most counted, A, C, G then T
        let mut best = A;
        for k in A + 1..=T {
            if counts[k] > counts[best] {
                best = k;
            }
        }
        let (mut count, mut symbol) = if counts[best] > 0 {
            (counts[best], 1usize << best)
        } else {
            (0, 0)
        };
        if count == 0 && counts[N] > 0 {
            (count, symbol) = (counts[N], 15);
        }
        if count >= counts[GAP] {
            consensus_row[i] = SYMBOLS[symbol];
            consensus.push(SYMBOLS[symbol]);
        } else {
            consensus_row[i] = b'-';
        }
    }
    Alignment {
        rows,
        consensus_row,
        consensus,
        profile,
    }
}
