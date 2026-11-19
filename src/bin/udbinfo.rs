use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use std::fs::File;
use std::io;
use std::io::{Read, Write};

use clap::{Arg, Command};

use rsearch::udb::Header;

/// prints an error and exits.
fn fail(message: &str) -> ! {
    eprintln!("Error: {message}");
    std::process::exit(1);
}

fn main() {
    let matches = Command::new("udbinfo")
        .about("Report a UDB file's header.")
        .long_about(
            "Report a UDB file's header.\n\n\
             Writes vsearch's --udbinfo report: the number of sequences, the word \
             length and the other fields of the header, as the file holds them. Like \
             vsearch, only the header is checked, and that the file is long enough for \
             its sequences.",
        )
        .arg(
            Arg::new("in")
                .short('i')
                .long("in")
                .value_name("INPUT")
                .help("Input UDB file, default to stdin"),
        )
        .arg(
            Arg::new("out")
                .short('o')
                .long("out")
                .value_name("OUTPUT")
                .help("Output report, default to stdout"),
        )
        .get_matches();

    let mut fin: Box<dyn Read> = match matches.get_one::<String>("in").map(String::as_str) {
        Some(f) => {
            Box::new(File::open(f).unwrap_or_else(|e| fail(&format!("cannot open {f}: {e}"))))
        }
        None => Box::new(io::stdin().lock()),
    };
    // the header, then the rest of the input counted, for its size
    let mut bytes = Vec::with_capacity(200);
    let mut size = fin
        .by_ref()
        .take(200)
        .read_to_end(&mut bytes)
        .unwrap_or_else(|e| fail(&format!("cannot read: {e}"))) as u64;
    size +=
        io::copy(&mut fin, &mut io::sink()).unwrap_or_else(|e| fail(&format!("cannot read: {e}")));
    let header = Header::parse(&bytes, size).unwrap_or_else(|e| fail(&e.to_string()));

    let mut fout: Box<dyn Write> = match matches.get_one::<String>("out").map(String::as_str) {
        Some(f) => {
            Box::new(File::create(f).unwrap_or_else(|e| fail(&format!("cannot create {f}: {e}"))))
        }
        None => Box::new(io::stdout()),
    };
    fout.write_all(header.report().as_bytes())
        .and_then(|()| fout.flush())
        .unwrap_or_else(|e| fail(&format!("cannot write: {e}")));
}
