use bio::io::fasta::Record;
#[derive(Debug, Default)]
/// A struct to store the size and annotation of sequences.
pub struct Cluster {
    pub size: u64,
    pub merged: Vec<String>,
}

impl Cluster {
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
        Cluster { size, merged }
    }
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
            1 => if name.is_none() && !pair[0].is_empty() {
                name = Some(pair[0]);
            },
            2 => if size.is_none() && pair[0] == "size" {
                if let Ok(x) = pair[1].parse().map(|x| if x > 0 { x } else { 1 }) {
                    size = Some(x);
                };
            },
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
}
