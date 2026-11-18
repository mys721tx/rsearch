use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read, Write};

use bio::io::{fasta, fastq};
use clap::{Arg, ArgAction, Command, value_parser};

use rsearch::cluster::{Cluster, Clusters, sort_by_size};

fn main() {
    let matches = Command::new("DeRep")
        .about("Remove duplications in FASTA or FASTQ sequences.")
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input FASTA or FASTQ file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file in the input format, default to stdout"),
        )
        .arg(
            Arg::new("strand")
                .short('s')
                .long("strand")
                .action(ArgAction::SetTrue)
                .help("Also merge reverse-complemented duplicates"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .value_name("MIN_SIZE")
                .value_parser(value_parser!(u64))
                .help("Minimum size of output sequences, disabled by default"),
        )
        .arg(
            Arg::new("max")
                .long("max")
                .value_name("MAX_SIZE")
                .value_parser(value_parser!(u64))
                .help("Maximum size of output sequences, disabled by default"),
        )
        .arg(
            Arg::new("top")
                .long("top")
                .value_name("TOP")
                .value_parser(value_parser!(usize))
                .help("Number of sequences to output, disabled by default"),
        )
        .arg(
            Arg::new("unsorted")
                .long("unsorted")
                .action(ArgAction::SetTrue)
                .help("Output in order of first appearance instead of decreasing size"),
        )
        .arg(
            Arg::new("relabel")
                .long("relabel")
                .value_name("PREFIX")
                .help("Relabel output sequences as PREFIX1, PREFIX2, ..."),
        )
        .arg(
            Arg::new("map")
                .long("map")
                .value_name("MAP")
                .help("Output file mapping each merged sequence to its output sequence"),
        )
        .get_matches();

    let fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => Box::new(File::open(f).unwrap()),
        None => Box::new(io::stdin()),
    };

    let mut fin = BufReader::new(fin);

    // FASTQ records start with '@', FASTA records with '>'.
    let is_fastq = fin.fill_buf().unwrap().first() == Some(&b'@');

    let strand = matches.get_flag("strand");

    let mut seqs = Clusters::default();

    if is_fastq {
        for record in fastq::Reader::from_bufread(fin).records() {
            let record = record.unwrap();
            if !record.seq().is_empty() {
                seqs.add(Cluster::build_from_fastq_record(&record), strand);
            }
        }
    } else {
        for record in fasta::Reader::from_bufread(fin).records() {
            let record = record.unwrap();
            if !record.seq().is_empty() {
                seqs.add(Cluster::build_from_record(&record), strand);
            }
        }
    }

    let mut clusters = seqs.into_clusters();

    if !matches.get_flag("unsorted") {
        sort_by_size(&mut clusters);
    }

    let min = matches.get_one::<u64>("min").copied().unwrap_or(0);
    let max = matches.get_one::<u64>("max").copied().unwrap_or(u64::MAX);
    let top = matches
        .get_one::<usize>("top")
        .copied()
        .unwrap_or(usize::MAX);

    let selected = clusters
        .iter()
        .filter(|c| c.size >= min && c.size <= max)
        .take(top);

    let relabel = matches.get_one::<String>("relabel");

    let label = |i: usize, c: &Cluster| match relabel {
        Some(prefix) => format!("{}{}", prefix, i + 1),
        None => c.merged[0].clone(),
    };

    let fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(File::create(f).unwrap()),
        None => Box::new(io::stdout()),
    };

    if is_fastq {
        let mut writer = fastq::Writer::new(fout);
        for (i, c) in selected.clone().enumerate() {
            let name = format!("{};size={}", label(i, c), c.size);
            let qual = c.qual.as_deref().unwrap_or_default();
            writer.write(name.as_str(), None, &c.seq, qual).unwrap();
        }
        writer.flush().unwrap();
    } else {
        let mut writer = fasta::Writer::new(fout);
        for (i, c) in selected.clone().enumerate() {
            let name = format!("{};size={}", label(i, c), c.size);
            writer.write(name.as_str(), None, &c.seq).unwrap();
        }
        writer.flush().unwrap();
    }

    if let Some(f) = matches.get_one::<String>("map") {
        let mut map = io::BufWriter::new(File::create(f).unwrap());
        for (i, c) in selected.enumerate() {
            let name = label(i, c);
            for member in &c.merged {
                writeln!(map, "{member}\t{name}").unwrap();
            }
        }
        map.flush().unwrap();
    }
}
