extern crate bio;
extern crate clap;

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{Read, Write};

use bio::io::fasta::{Reader, Record, Writer};
use clap::{App, Arg};

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
    let matches = App::new("DeRep")
        .about("Remove duplications in FASTA sequences.")
        .arg(
            Arg::with_name("in")
                .short("i")
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA file, default to stdin")
                .takes_value(true),
        ).arg(
            Arg::with_name("out")
                .short("o")
                .long("out")
                .value_name("OUTPUT")
                .help("Output FASTA file, default to stdout")
                .takes_value(true),
        ).get_matches();

    let fin: Box<Read> = match matches.value_of("in") {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let reader = Reader::new(fin);

    let mut seqs: HashMap<Vec<u8>, Cluster> = HashMap::new();

    for record in reader.records() {
        let record = record.unwrap();

        let key = record.seq().to_vec();

        let cluster = Cluster::build_from_record(&record);

        let entry = seqs.entry(key).or_insert_with(Cluster::new);

        entry.size += cluster.size;
        entry.merged.extend(cluster.merged);
    }

    let fout: Box<Write> = match matches.value_of("out") {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    let mut writer = Writer::new(fout);

    for (k, v) in &seqs {
        let name = format!("{};size={}", v.merged[0], v.size);
        writer.write(name.as_str(), Option::None, k).unwrap();
    }
}
