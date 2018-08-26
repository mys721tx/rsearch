extern crate bio;

use bio::io::fasta::{Reader, Writer, Record};

use std::collections::HashMap;
use std::io;

#[derive(Debug)]
struct Cluster {
    size: u64,
    merged: Vec<String>,
}

impl Cluster {
    fn new() -> Cluster {
        let merged: Vec<String> = Vec::new();
        Cluster { size: 0, merged }
    }
    fn build_from_record(record: &Record) -> Cluster {
        let (name, size) = parse_id(record.id());

        let merged = vec![name.unwrap_or_default("sequence").to_string()];
        Cluster {
            size: size.unwrap_or_default(1),
            merged,
        }
    }
}

enum Field<T> {
    Set(T),
    Unset,
}

impl<T> Field<T> {
    fn unwrap_or_default(self, default: T) -> T {
        match self {
            Field::Set(x) => x,
            Field::Unset => default,
        }
    }
}

fn parse_id(id: &str) -> (Field<&str>, Field<u64>) {
    let mut name = Field::Unset;
    let mut size: Field<u64> = Field::Unset;

    for item in id.split(';') {
        let pair: Vec<&str> = item.split('=').collect();
        match pair.len() {
            1 => match name {
                Field::Unset => name = Field::Set(pair[0]),
                Field::Set(_) => (),
            },
            2 => match size {
                Field::Unset => if pair[0] == "size" {
                    if let Ok(x) = pair[1].parse() {
                        if x > 0 {
                            size = Field::Set(x)
                        }
                    }
                },
                Field::Set(_) => (),
            },
            _ => (),
        }
    }

    (name, size)
}

fn main() {
    let mut seqs: HashMap<Vec<u8>, Cluster> = HashMap::new();

    let reader = Reader::new(io::stdin());

    for record in reader.records() {
        let record = record.unwrap();

        let key = record.seq().to_vec();

        let cluster = Cluster::build_from_record(&record);

        let entry = seqs.entry(key).or_insert_with(Cluster::new);

        entry.size += cluster.size;
        entry.merged.extend(cluster.merged);
    }

    let mut writer = Writer::new(io::stdout());

    for (k, v) in &seqs {
        let name = format!("{};size={}", v.merged[0], v.size);
        writer.write(name.as_str(), Option::None, k).unwrap();
    }
}
