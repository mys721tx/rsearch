extern crate bio;

use bio::io::fasta;

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
    fn build_from_record(record: &bio::io::fasta::Record) -> Cluster {
        let merged = vec![record.id().to_string()];

        Cluster { size: 1, merged }
    }
}

fn main() {
    let mut seqs: HashMap<Vec<u8>, Cluster> = HashMap::new();

    let reader = fasta::Reader::new(io::stdin());

    for record in reader.records() {
        let record = record.unwrap();

        let key = record.seq().to_vec();

        let cluster = Cluster::build_from_record(&record);

        let entry = seqs.entry(key).or_insert_with(Cluster::new);

        entry.size += cluster.size;
        entry.merged.extend(cluster.merged);
    }

    println!("{:?}", seqs);
}
