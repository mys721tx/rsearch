use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Write};
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{Arg, ArgAction, ArgGroup, Command};

use rsearch::hits::UcQueries;
use rsearch::label::abundance;
use rsearch::otutable::OtuTable;

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

/// the time now, as `YYYY-MM-DDTHH:MM:SS` in UTC.
fn now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs()) as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    // days to a civil date (Howard Hinnant's algorithm)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

fn main() {
    let matches = Command::new("otutab")
        .about("Count hits per sample and OTU: classic, mothur or BIOM table.")
        .long_about(
            "Count hits per sample and OTU: classic, mothur or BIOM table.\n\n\
             As vsearch's --otutabout, --mothur_shared_out and --biomout: each query of \
             the hits counts its abundance (;size=) for its sample (sample= or \
             barcodelabel=, or its label's leading letters, digits and underscores) and \
             the OTU of its first hit's target (otu=, or its label up to the first \
             semicolon). Every database sequence is an OTU, with its tax= as taxonomy.",
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
                .help("Output table, default to stdout"),
        )
        .arg(
            Arg::new("mothur")
                .long("mothur")
                .action(ArgAction::SetTrue)
                .help("Write a mothur shared table"),
        )
        .arg(
            Arg::new("biom")
                .long("biom")
                .action(ArgAction::SetTrue)
                .help("Write a BIOM 1.0 table"),
        )
        .group(ArgGroup::new("format").args(["mothur", "biom"]))
        .get_matches();

    let db: Vec<String> = rsearch::database::read(matches.get_one::<String>("db").unwrap())
        .unwrap_or_else(|e| fail(&e))
        .into_iter()
        .map(|(header, _)| header)
        .collect();
    let hits_path = matches.get_one::<String>("hits").unwrap();
    let hits_file =
        File::open(hits_path).unwrap_or_else(|e| fail(&format!("cannot open {hits_path}: {e}")));
    let mut uc = UcQueries::new(BufReader::new(hits_file));
    let mut table = OtuTable::default();
    let mut matched = vec![false; db.len()];
    let target = |t: u64| db.get(t as usize).map(String::as_str);
    while let Some(query) = uc.next_query(target) {
        let (label, hits) = query.unwrap_or_else(|e| fail(&e));
        let size = abundance(&label);
        table.add(
            Some(&label),
            hits.first().map(|h| db[h.target as usize].as_str()),
            size,
        );
        for hit in &hits {
            matched[hit.target as usize] = true;
        }
    }
    // the database sequences without hits, as vsearch adds them last
    for (label, _) in db.iter().zip(&matched).filter(|(_, m)| !**m) {
        table.add(None, Some(label), 0);
    }

    let out_path = matches.get_one::<String>("out");
    let text = if matches.get_flag("mothur") {
        table.mothur()
    } else if matches.get_flag("biom") {
        table.biom(
            out_path.map_or("-", String::as_str),
            &format!("rsearch {}", env!("CARGO_PKG_VERSION")),
            &now(),
        )
    } else {
        table.classic()
    };
    let mut out: Box<dyn Write> = match out_path {
        Some(f) => Box::new(BufWriter::new(
            File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))),
        )),
        None => Box::new(BufWriter::new(io::stdout())),
    };
    out.write_all(text.as_bytes())
        .and_then(|()| out.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
