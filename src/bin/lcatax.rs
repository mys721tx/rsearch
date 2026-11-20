use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Write};

use clap::{Arg, Command};

use rsearch::hits::UcQueries;
use rsearch::taxonomy::lca;

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// a share above 0.5, at most 1.
fn cutoff(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(x) if x > 0.5 && x <= 1.0 => Ok(x),
        _ => Err("expected a number above 0.5, at most 1".to_string()),
    }
}

fn main() {
    let matches = Command::new("lcatax")
        .about("Last common ancestor of each query's top hits.")
        .long_about(
            "Last common ancestor of each query's top hits.\n\n\
             As vsearch's --lcaout: for each query of the hits, its label, a tab, and \
             the ;tax= levels of the database sequences it hits agree on, level by level \
             while at least --cutoff of them agree, as vsearch's majority vote picks the \
             name.",
        )
        .arg(
            Arg::new("hits")
                .value_name("HITS")
                .required(true)
                .help("UC file of hits"),
        )
        .arg(
            Arg::new("db")
                .long("db")
                .value_name("FILE")
                .required(true)
                .help("Database the hits were found in: FASTA, FASTQ or UDB file"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output file, default to stdout"),
        )
        .arg(
            Arg::new("cutoff")
                .long("cutoff")
                .value_name("SHARE")
                .value_parser(cutoff)
                .default_value("1")
                .help("Share of the hits that must agree on a level, above 0.5"),
        )
        .get_matches();
    let cutoff = *matches.get_one::<f64>("cutoff").unwrap();

    let db: Vec<String> = rsearch::database::read(matches.get_one::<String>("db").unwrap())
        .unwrap_or_else(|e| fail(&e))
        .into_iter()
        .map(|(header, _)| header)
        .collect();
    let hits_path = matches.get_one::<String>("hits").unwrap();
    let hits_file =
        File::open(hits_path).unwrap_or_else(|e| fail(&format!("cannot open {hits_path}: {e}")));
    let mut uc = UcQueries::new(BufReader::new(hits_file));
    let mut out: Box<dyn Write> = match matches.get_one::<String>("out") {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let target = |t: u64| db.get(t as usize).map(String::as_str);
    while let Some(query) = uc.next_query(target) {
        let (label, hits) = query.unwrap_or_else(|e| fail(&e));
        let headers: Vec<&[u8]> = hits
            .iter()
            .map(|h| db[h.target as usize].as_bytes())
            .collect();
        out.write_all(label.as_bytes())
            .and_then(|()| out.write_all(b"\t"))
            .and_then(|()| out.write_all(&lca(&headers, cutoff)))
            .and_then(|()| out.write_all(b"\n"))
            .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
    }
    out.flush()
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
