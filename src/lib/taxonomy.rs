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
