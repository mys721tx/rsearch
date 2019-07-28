use jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{Read, Write};

use bio::io::fasta::{Reader, Writer};
use clap::{App, Arg};

use rsearch::cluster::Cluster;

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

    let fin: Box<dyn Read> = match matches.value_of("in") {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let reader = Reader::new(fin);

    let mut seqs: HashMap<Vec<u8>, Cluster> = HashMap::new();

    for record in reader.records() {
        let record = record.unwrap();

        let key = record.seq().to_vec();

        let cluster = Cluster::build_from_record(&record);

        let entry = seqs.entry(key).or_default();

        entry.size += cluster.size;
        entry.merged.extend(cluster.merged);
    }

    let fout: Box<dyn Write> = match matches.value_of("out") {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    let mut writer = Writer::new(fout);

    for (k, v) in &seqs {
        let name = format!("{};size={}", v.merged[0], v.size);
        writer.write(name.as_str(), Option::None, k).unwrap();
    }
}
