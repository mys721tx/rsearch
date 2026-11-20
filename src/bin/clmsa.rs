use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use clap::{Arg, ArgAction, ArgGroup, Command};

use rsearch::cluster::reverse_complement_keeping_case;
use rsearch::fastx::{Reader, Record};
use rsearch::hits::ClusterRecord;
use rsearch::label::{Label, checked_abundance};
use rsearch::msa::{Member, align, member_ops};

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// A record of a cluster, with its alignment to the centroid.
struct Entry {
    header: String,
    seq: Vec<u8>,
    size: u64,
    minus: bool,
    ops: Vec<(u8, usize)>,
}

fn main() {
    let matches = Command::new("clmsa")
        .about("Align each cluster's members to its centroid, from clustering hits.")
        .long_about(
            "Align each cluster's members to its centroid, from clustering hits.\n\n\
             As vsearch --msaout, --consout and --profile: the S and H records of the \
             hits of clust or unoise give the cluster of each record of the input \
             clustered, in order, and each member's alignment to its centroid. Each \
             cluster's records are aligned together, members of the minus strand \
             reverse-complemented, and its profile counts each column's residues, \
             weighted by abundance (;size=). Writes one of the multiple alignments with \
             their consensus (--msa, the default), the consensus sequences \
             (--consensus) or the profiles (--profile), FASTA unwrapped.",
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
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file, default to stdout"),
        )
        .arg(
            Arg::new("msa")
                .long("msa")
                .action(ArgAction::SetTrue)
                .help("Write each cluster's multiple alignment and consensus, as vsearch --msaout"),
        )
        .arg(
            Arg::new("consensus")
                .long("consensus")
                .action(ArgAction::SetTrue)
                .help("Write each cluster's consensus sequence, as vsearch --consout"),
        )
        .arg(
            Arg::new("profile")
                .long("profile")
                .action(ArgAction::SetTrue)
                .help("Write each cluster's profile, as vsearch --profile"),
        )
        .group(ArgGroup::new("format").args(["msa", "consensus", "profile"]))
        .arg(
            Arg::new("no-sizein")
                .long("no-sizein")
                .action(ArgAction::SetTrue)
                .help("Count each sequence once, rather than by its abundance (;size=)"),
        )
        .arg(
            Arg::new("no-sizeout")
                .long("no-sizeout")
                .action(ArgAction::SetTrue)
                .help("Do not annotate the labels with abundances (;size=)"),
        )
        .get_matches();
    let sizein = !matches.get_flag("no-sizein");
    let sizeout = !matches.get_flag("no-sizeout");

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

    // every record is matched with its cluster and checked before anything
    // is written
    let mut clusters: Vec<Vec<Entry>> = Vec::new();
    for line in hits.lines() {
        let line = line.unwrap_or_else(|e| fail(&format!("cannot read {path}: {e}")));
        let Some(hit) = ClusterRecord::parse(&line, clusters.len()).unwrap_or_else(|e| fail(&e))
        else {
            continue;
        };
        if !next_record(&mut record) {
            fail(&format!("no record for the hit of {}", hit.label));
        }
        let header = record.header();
        if header != hit.label {
            fail(&format!(
                "the record {header} does not match the hit of {}",
                hit.label
            ));
        }
        let size =
            checked_abundance(header).unwrap_or_else(|e| fail(&format!("{e} in read {header}")));
        let ops = if hit.cluster == clusters.len() {
            clusters.push(Vec::new());
            Vec::new()
        } else {
            let centroid = clusters[hit.cluster][0].seq.len();
            member_ops(hit.cigar, record.seq().len(), centroid)
                .unwrap_or_else(|e| fail(&format!("{e}, in the hit of {header}")))
        };
        clusters[hit.cluster].push(Entry {
            header: header.to_string(),
            seq: record.seq().to_vec(),
            size,
            minus: hit.minus,
            ops,
        });
    }
    if next_record(&mut record) {
        fail(&format!("no hit for the record {}", record.header()));
    }

    let mut out: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let label = Label {
        sizeout,
        ..Label::default()
    };
    let mut text = String::new();
    for entries in &clusters {
        let weight = |e: &Entry| if sizein { e.size } else { 1 };
        let total: u64 = entries.iter().map(weight).sum();
        let members: Vec<Member> = entries[1..]
            .iter()
            .map(|e| {
                let mut seq = Vec::new();
                if e.minus {
                    reverse_complement_keeping_case(&mut seq, &e.seq);
                } else {
                    seq.clone_from(&e.seq);
                }
                Member {
                    seq,
                    weight: weight(e),
                    ops: e.ops.clone(),
                }
            })
            .collect();
        let centroid = &entries[0];
        let aligned = align(&centroid.seq, weight(centroid), &members);
        let cluster_label = Label {
            cluster_size: Some(entries.len() as u64),
            ..label.clone()
        }
        .format(&centroid.header, &aligned.consensus, total, 0, None);
        text.clear();
        if matches.get_flag("consensus") {
            text.push_str(&format!(
                ">centroid={cluster_label}\n{}\n",
                String::from_utf8_lossy(&aligned.consensus)
            ));
        } else if matches.get_flag("profile") {
            // an empty sequence line, as vsearch with --fasta_width 0
            text.push_str(&format!(">centroid={cluster_label}\n\n"));
            for i in 0..aligned.consensus_row.len() {
                text.push_str(&aligned.profile_line(i));
                text.push('\n');
            }
            text.push('\n');
        } else {
            text.push('\n');
            for (i, (e, row)) in entries.iter().zip(&aligned.rows).enumerate() {
                let prefix = if i == 0 { "*" } else { "" };
                let header = label.format(&e.header, &e.seq, e.size, 0, None);
                text.push_str(&format!(
                    ">{prefix}{header}\n{}\n",
                    String::from_utf8_lossy(row)
                ));
            }
            text.push_str(&format!(
                ">consensus\n{}\n",
                String::from_utf8_lossy(&aligned.consensus_row)
            ));
        }
        out.write_all(text.as_bytes())
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    out.flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
