use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use clap::parser::ValueSource;
use clap::{Arg, Command};

use rsearch::clustering::{cluster_args, run};

/// a finite number.
fn finite(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(x) if x.is_finite() => Ok(x),
        _ => Err("expected a number".to_string()),
    }
}

fn main() {
    let command = Command::new("unoise")
        .about("Denoise amplicons with the UNOISE3 rule; write centroids and hits.")
        .long_about(
            "Denoise amplicons with the UNOISE3 rule; write centroids and hits.\n\n\
             As vsearch --cluster_unoise on input sorted by abundance (sortby --size), of \
             abundance at least --minsize (sqfilt --min-size): each sequence is searched \
             against the centroids so far, as by clust, and joins the cluster of its best \
             hit of at least --id identity whose target is at least 2^(alpha d + 1) times \
             as abundant, for d mismatches (or that has none), or becomes a centroid. The \
             centroids and hits are written as by clust.",
        )
        .arg(
            Arg::new("alpha")
                .long("alpha")
                .value_name("ALPHA")
                .value_parser(finite)
                .default_value("2")
                .allow_negative_numbers(true)
                .help("UNOISE's alpha: how much more abundant a target is, a mismatch"),
        );
    let mut command = cluster_args(command).mut_arg("id", |a| {
        a.default_value("0.9")
            .help("Accept hits of at least this identity, 0 to 1")
    });
    let matches = command.clone().get_matches();
    // the identity is UNOISE's minimum, vsearch's --weak_id for it
    if matches.value_source("weak-id") == Some(ValueSource::CommandLine) {
        command
            .error(
                clap::error::ErrorKind::ArgumentConflict,
                "unoise takes its minimum identity as --id, not --weak-id",
            )
            .exit();
    }
    let alpha = *matches.get_one::<f64>("alpha").unwrap();
    run(&matches, Some(alpha));
}
