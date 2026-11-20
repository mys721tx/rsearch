use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read};

use clap::{Arg, Command};

use rsearch::fastx::{Reader, Record};
use rsearch::writer::Writer;

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// A record of a cluster: header, sequence and scores.
struct Member {
    header: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// the cluster of an S or H record of `line` and its query label, or
/// `None` for a C record.
fn cluster_record(line: &str, clusters: usize) -> Result<Option<(usize, &str)>, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() != 10 {
        return Err(format!("not a UC record: {line}"));
    }
    let kind = fields[0];
    if kind == "C" {
        return Ok(None);
    }
    let cluster: usize = fields[1]
        .parse()
        .map_err(|_| format!("invalid cluster number in UC record: {line}"))?;
    match kind {
        // a centroid opens the next cluster
        "S" if cluster == clusters => Ok(Some((cluster, fields[8]))),
        "H" if cluster < clusters => Ok(Some((cluster, fields[8]))),
        "S" | "H" => Err(format!(
            "cluster {cluster} out of order in UC record: {line}"
        )),
        _ => Err(format!("not a cluster record: {line}")),
    }
}

fn main() {
    let matches = Command::new("clspl")
        .about("Split records into one file per cluster, from clustering hits.")
        .long_about(
            "Split records into one file per cluster, from clustering hits.\n\n\
             As vsearch --clusters: the S and H records of the hits of clust or unoise \
             give the cluster of each record of the input clustered, in order; each \
             cluster's records are written, in input order and in the input format, \
             to PREFIX followed by the cluster's number, from 0.",
        )
        .arg(
            Arg::new("hits")
                .long("hits")
                .value_name("HITS")
                .required(true)
                .help("UC file of the clusters"),
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA or FASTQ file clustered, default to stdin"),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .value_name("PREFIX")
                .required(true)
                .help("Prefix of the files written, followed by each cluster's number"),
        )
        .get_matches();

    let path = matches.get_one::<String>("hits").unwrap();
    let hits = BufReader::new(
        File::open(path).unwrap_or_else(|e| fail(&format!("cannot open {path}: {e}"))),
    );
    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin()),
    };
    let mut reader = Reader::new(BufReader::new(fin));
    let mut record = Record::new();
    let mut next_record = |record: &mut Record| {
        reader.read(record).unwrap_or_else(|e| fail(&e.to_string()));
        !record.is_empty()
    };

    // every record is matched with its cluster before any file is written
    let mut clusters: Vec<Vec<Member>> = Vec::new();
    let mut is_fastq = false;
    for line in hits.lines() {
        let line = line.unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
        let Some((cluster, label)) =
            cluster_record(&line, clusters.len()).unwrap_or_else(|e| fail(&e))
        else {
            continue;
        };
        if !next_record(&mut record) {
            fail(&format!("no record for the hit of {label}"));
        }
        if record.header() != label {
            fail(&format!(
                "the record {} does not match the hit of {label}",
                record.header()
            ));
        }
        is_fastq = record.is_fastq();
        if cluster == clusters.len() {
            clusters.push(Vec::new());
        }
        clusters[cluster].push(Member {
            header: record.header().to_string(),
            seq: record.seq().to_vec(),
            qual: record.qual().to_vec(),
        });
    }
    if next_record(&mut record) {
        fail(&format!("no hit for the record {}", record.header()));
    }

    let prefix = matches.get_one::<String>("prefix").unwrap();
    for (i, members) in clusters.iter().enumerate() {
        let name = format!("{prefix}{i}");
        let file =
            File::create(&name).unwrap_or_else(|e| fail(&format!("cannot create {name}: {e}")));
        let mut writer = Writer::new(BufWriter::new(file), !is_fastq, 0);
        for m in members {
            writer
                .write(&m.header, None, &m.seq, &m.qual)
                .unwrap_or_else(|e| fail(&format!("cannot write {name}: {e}")));
        }
        writer
            .flush()
            .unwrap_or_else(|e| fail(&format!("cannot write {name}: {e}")));
    }
}
