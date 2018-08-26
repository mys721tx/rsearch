use bio::io::fasta::Record;
#[derive(Debug, Default)]
pub struct Cluster {
    pub size: u64,
    pub merged: Vec<String>,
}

impl Cluster {
    pub fn new() -> Cluster {
        let merged: Vec<String> = Vec::new();
        Cluster { size: 0, merged }
    }
    pub fn build_from_record(record: &Record) -> Cluster {
        let (name, size) = parse_id(record.id());

        let merged = vec![name.to_string()];
        Cluster { size, merged }
    }
}

pub fn parse_id(id: &str) -> (&str, u64) {
    let mut name = None;
    let mut size: Option<u64> = None;

    for item in id.split(';') {
        let pair: Vec<&str> = item.split('=').collect();
        match pair.len() {
            1 => if name.is_none() {
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
        assert_eq!(
            parse_id(input),
            expected,
            "parse_id should parse a semicolon seperated list."
        );

        input = "";
        expected = ("sequence", 1);
        assert_eq!(
            parse_id(input),
            expected,
            "parse_id should return name of \"sequence\" and size of 1 by default."
        );

        input = "egg;spam=fried;";
        expected = ("egg", 1);
        assert_eq!(
            parse_id(input),
            expected,
            "parse_id should return a size of 1 if not such field is in the input."
        );

        input = "species=9606;spam=fried;size=100";
        expected = ("sequence", 100);
        assert_eq!(
            parse_id(input),
            expected,
            "parse_id should return a name of \"sequence\" of 1 if there is no non-empty field with an equal sign."
        );

        input = "maple;lily;iris";
        expected = ("maple", 1);
        assert_eq!(
            parse_id(input),
            expected,
            "parse_id should return the first non-empty field with an equal sign as the name."
        );

        input = "size=1;size=2;size=3;size=4";
        expected = ("sequence", 1);
        assert_eq!(
            parse_id(input),
            expected,
            "parse_id should return the value of the first non-empty field with an equal sign as size."
        );

        input = "size=abc;size=def;size=ghi";
        expected = ("sequence", 1);
        assert_eq!(
            parse_id(input),
            expected,
            "parse_id should return a size of 1 when all such fields are not numbers."
        );
    }
}
