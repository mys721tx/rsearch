//! Taxonomy annotations: the `;tax=` of a header, as vsearch parses it,
//! split into its nine levels, `d:domain,k:kingdom,p:phylum,c:class,
//! o:order,f:family,g:genus,s:species,t:strain`.

use crate::label::find_attribute_with;

/// The level letters, in order.
pub const LEVELS: [u8; 9] = *b"dkpcofgst";

/// the name at each level of the first `tax=` annotation of `header`, as
/// byte ranges of the header, empty for the levels it lacks. A level
/// letter, in either case, and a colon start a name, which runs to the
/// next comma; anything else up to a comma is skipped, and a later name
/// of a level replaces an earlier one.
///
/// # Example
/// ```rust
/// # use rsearch::taxonomy::split;
/// let header = b"r1;tax=d:Bacteria,G:Bacillus,x:y,g;size=2";
/// let levels = split(header);
/// assert_eq!(&header[levels[0].clone()], b"Bacteria");
/// assert_eq!(&header[levels[6].clone()], b"Bacillus");
/// assert!(levels[2].is_empty());
/// ```
pub fn split(header: &[u8]) -> [std::ops::Range<usize>; 9] {
    let mut levels: [std::ops::Range<usize>; 9] = Default::default();
    let Some((start, end)) = find_attribute_with(header, b"tax=", |c| c != b';', true) else {
        return levels;
    };
    let mut offset = start + 4;
    while offset < end {
        let comma = header[offset..end]
            .iter()
            .position(|&c| c == b',')
            .map_or(end, |i| offset + i);
        if let Some(level) = LEVELS
            .iter()
            .position(|&l| l == header[offset].to_ascii_lowercase())
            && offset + 1 < header.len()
            && header[offset + 1] == b':'
        {
            // the name runs from after the colon to the comma, or the end
            let name_end = if offset + 2 <= end {
                header[offset + 2..end]
                    .iter()
                    .position(|&c| c == b',')
                    .map_or(end, |i| offset + 2 + i)
            } else {
                end
            };
            levels[level] = offset + 2..name_end.max(offset + 2);
        }
        offset = if comma < end { comma + 1 } else { end };
    }
    levels
}

/// the last common ancestor of the taxonomies in `headers`, as vsearch's
/// `--lcaout`. At each level, a candidate is voted for by the headers
/// agreeing with it on that level and those above (a majority vote, as
/// vsearch's). The levels are written while at least `cutoff` of the
/// headers agree with their candidate, empty names skipped.
///
/// # Example
/// ```rust
/// # use rsearch::taxonomy::lca;
/// let headers: [&[u8]; 3] = [b"a;tax=d:B,p:P1,g:G1", b"b;tax=d:B,p:P1,g:G2", b"c;tax=d:B,p:P2"];
/// assert_eq!(lca(&headers, 0.6), b"d:B,p:P1");
/// assert_eq!(lca(&headers, 1.0), b"d:B");
/// ```
pub fn lca(headers: &[&[u8]], cutoff: f64) -> Vec<u8> {
    let names: Vec<[&[u8]; 9]> = headers
        .iter()
        .map(|header| split(header).map(|range| &header[range]))
        .collect();
    // the candidate of each level, and its votes
    let mut candidate = [0usize; 9];
    let mut votes = [0i64; 9];
    let agree = |a: usize, b: usize, k: usize| (0..=k).all(|j| names[a][j] == names[b][j]);
    for i in 0..names.len() {
        for k in 0..LEVELS.len() {
            if votes[k] == 0 {
                candidate[k] = i;
                votes[k] = 1;
            } else if agree(candidate[k], i, k) {
                votes[k] += 1;
            } else {
                votes[k] -= 1;
            }
        }
    }
    let mut out = Vec::new();
    if names.is_empty() {
        return out;
    }
    for k in 0..LEVELS.len() {
        let matching = (0..names.len())
            .filter(|&i| agree(candidate[k], i, k))
            .count();
        if (matching as f64) / (names.len() as f64) < cutoff {
            break;
        }
        let name = names[candidate[k]][k];
        if !name.is_empty() {
            if !out.is_empty() {
                out.push(b',');
            }
            out.push(LEVELS[k]);
            out.push(b':');
            out.extend_from_slice(name);
        }
    }
    out
}
