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
