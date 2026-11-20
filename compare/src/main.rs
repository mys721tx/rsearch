//! `rscompare`: compares rsearch with vsearch. See `compare/README.md`.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Arg, ArgAction, Command};
use rscompare::{runner, sources};

fn cli() -> Command {
    Command::new("rscompare")
        .about("Compares rsearch with vsearch")
        .subcommand_required(true)
        .subcommand(
            Command::new("run")
                .about(
                    "Runs the sources' cases and prints a DIFF line for each failure and the tally",
                )
                .arg(
                    Arg::new("sources")
                        .num_args(1..)
                        .required(true)
                        .help("Sources, or all"),
                )
                .arg(
                    Arg::new("seeds")
                        .long("seeds")
                        .help("Seeds such as 1-50,70 (default: each source's own)"),
                )
                .arg(
                    Arg::new("jobs")
                        .long("jobs")
                        .short('j')
                        .value_parser(clap::value_parser!(usize))
                        .help("Cases run at once (default: the cores)"),
                )
                .arg(
                    Arg::new("vsearch")
                        .long("vsearch")
                        .help("vsearch (default: compare/vsearch-build/vsearch)"),
                )
                .arg(
                    Arg::new("bin")
                        .long("bin")
                        .help("rsearch's tools (default: target/release)"),
                )
                .arg(
                    Arg::new("out")
                        .long("out")
                        .help("Append the records, a line of JSON each, to this file"),
                )
                .arg(
                    Arg::new("tmp")
                        .long("tmp")
                        .help("Directory for the cases' files (default: $TMPDIR/rscompare)"),
                )
                .arg(
                    Arg::new("keep")
                        .long("keep")
                        .action(ArgAction::SetTrue)
                        .help("Keep the files of passing cases too"),
                ),
        )
        .subcommand(
            Command::new("time")
                .about("Times two builds of the tools on the benchmark inputs, alternating runs")
                .arg(
                    Arg::new("baseline")
                        .required(true)
                        .help("The baseline build's tools"),
                )
                .arg(Arg::new("new").required(true).help("The new build's tools"))
                .arg(
                    Arg::new("runs")
                        .long("runs")
                        .default_value("15")
                        .value_parser(clap::value_parser!(usize)),
                )
                .arg(
                    Arg::new("inputs")
                        .long("inputs")
                        .help("The benchmark inputs (default: target/tmp)"),
                ),
        )
        .subcommand(
            Command::new("build-vsearch")
                .about("Builds vsearch from its source tree, for the platform it runs on")
                .arg(
                    Arg::new("source")
                        .help("vsearch's repository (default: ../vsearch next to rsearch)"),
                )
                .arg(Arg::new("out").help("Output directory (default: compare/vsearch-build)"))
                .arg(
                    Arg::new("arch")
                        .long("arch")
                        .help("Architecture, as uname -m (default: this machine's)"),
                )
                .arg(
                    Arg::new("os")
                        .long("os")
                        .help("Operating system, as uname -s (default: this machine's)"),
                )
                .arg(
                    Arg::new("dry-run")
                        .long("dry-run")
                        .action(ArgAction::SetTrue)
                        .help("Print the chosen sources and flags, and build nothing"),
                ),
        )
        .subcommand(
            Command::new("fetch")
                .about(
                    "Downloads the real datasets of datasets.json and prepares their working files",
                )
                .arg(
                    Arg::new("ids")
                        .num_args(0..)
                        .help("Dataset ids (default: all)"),
                )
                .arg(
                    Arg::new("dir")
                        .long("dir")
                        .help("Data directory (default: $RSCOMPARE_DATA or $TMPDIR/rsearch-data)"),
                )
                .arg(
                    Arg::new("manifest")
                        .long("manifest")
                        .help("The manifest (default: compare/datasets.json)"),
                )
                .arg(
                    Arg::new("max-reads")
                        .long("max-reads")
                        .value_name("N")
                        .value_parser(clap::value_parser!(u64).range(1..))
                        .help("Keep the first N reads of each FASTQ file"),
                )
                .arg(
                    Arg::new("vsearch")
                        .long("vsearch")
                        .value_name("PATH")
                        .help("vsearch, for the steps that need it (UDB files)"),
                )
                .arg(flag(
                    "list",
                    "List the datasets and their sizes, fetch nothing",
                ))
                .arg(flag(
                    "download-only",
                    "Download and check, but do not unpack",
                ))
                .arg(flag(
                    "no-references",
                    "Do not add the reference databases the chosen datasets use",
                ))
                .arg(flag("force", "Rebuild the working files"))
                .arg(
                    Arg::new("retries")
                        .long("retries")
                        .default_value("5")
                        .value_parser(clap::value_parser!(u32)),
                )
                .arg(
                    Arg::new("quiet")
                        .short('q')
                        .long("quiet")
                        .action(ArgAction::SetTrue),
                ),
        )
        .subcommand(
            Command::new("list")
                .about("Lists the sources, their default seeds and what they cover"),
        )
        .subcommand(
            Command::new("summary")
                .about("Outcomes, coverage, times and caveats of recorded runs")
                .arg(Arg::new("results").num_args(1..).required(true)),
        )
        .subcommand(
            Command::new("translate")
                .about("Prints the rsearch pipeline for a vsearch command line")
                .arg(
                    Arg::new("workdir")
                        .long("workdir")
                        .default_value(".")
                        .help("Directory for the pipeline's files"),
                )
                .arg(
                    Arg::new("json")
                        .long("json")
                        .action(ArgAction::SetTrue)
                        .help("Print the pipeline as JSON"),
                )
                .arg(
                    Arg::new("vsearch")
                        .num_args(1..)
                        .trailing_var_arg(true)
                        .allow_hyphen_values(true)
                        .required(true)
                        .help("vsearch's arguments, after --"),
                ),
        )
}

fn flag(name: &'static str, help: &'static str) -> Arg {
    Arg::new(name)
        .long(name)
        .action(ArgAction::SetTrue)
        .help(help)
}

fn main() -> ExitCode {
    let m = cli().get_matches();
    match m.subcommand() {
        Some(("translate", m)) => translate(m),
        Some(("run", m)) => run(m),
        Some(("list", _)) => list(),
        Some(("time", m)) => {
            let inputs = m
                .get_one::<String>("inputs")
                .map_or_else(|| runner::repo_root().join("target/tmp"), PathBuf::from);
            let (base, new) = (
                m.get_one::<String>("baseline").unwrap(),
                m.get_one::<String>("new").unwrap(),
            );
            report(rscompare::timing::run(
                base.as_ref(),
                new.as_ref(),
                &inputs,
                *m.get_one::<usize>("runs").unwrap(),
            ))
        }
        Some(("build-vsearch", m)) => {
            let repo = runner::repo_root();
            let source = m
                .get_one::<String>("source")
                .map_or_else(|| repo.parent().unwrap().join("vsearch"), PathBuf::from);
            let out = m
                .get_one::<String>("out")
                .map_or_else(|| repo.join("compare/vsearch-build"), PathBuf::from);
            let mut opts = rscompare::vsearch_build::Options::default();
            if let Some(a) = m.get_one::<String>("arch") {
                opts.arch = a.clone();
            }
            if let Some(o) = m.get_one::<String>("os") {
                opts.os = o.clone();
            }
            opts.dry_run = m.get_flag("dry-run");
            match rscompare::vsearch_build::build(&source, &out, &opts) {
                Ok(bin) if !opts.dry_run => {
                    let v = std::process::Command::new(&bin).arg("--version").output();
                    if let Ok(v) = v {
                        print!(
                            "{}",
                            String::from_utf8_lossy(&v.stderr)
                                .lines()
                                .next()
                                .unwrap_or("")
                        );
                        println!();
                    }
                    ExitCode::SUCCESS
                }
                other => report(other.map(|_| ())),
            }
        }
        Some(("fetch", m)) => {
            use rscompare::fetch;
            let opts = fetch::Options {
                ids: m
                    .get_many::<String>("ids")
                    .map(|v| v.cloned().collect())
                    .unwrap_or_default(),
                dir: m
                    .get_one::<String>("dir")
                    .map_or_else(fetch::default_dir, PathBuf::from),
                manifest: m
                    .get_one::<String>("manifest")
                    .map_or_else(fetch::default_manifest, PathBuf::from),
                max_reads: m.get_one::<u64>("max-reads").copied(),
                vsearch: m.get_one::<String>("vsearch").map(PathBuf::from),
                list: m.get_flag("list"),
                download_only: m.get_flag("download-only"),
                no_references: m.get_flag("no-references"),
                force: m.get_flag("force"),
                retries: *m.get_one::<u32>("retries").unwrap(),
                quiet: m.get_flag("quiet"),
            };
            report(fetch::run(&opts))
        }
        Some(("summary", m)) => {
            let paths: Vec<PathBuf> = m
                .get_many::<String>("results")
                .unwrap()
                .map(PathBuf::from)
                .collect();
            match runner::summary(&paths) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!("rscompare: {e}");
                    ExitCode::FAILURE
                }
            }
        }
        _ => unreachable!(),
    }
}

fn translate(m: &clap::ArgMatches) -> ExitCode {
    let workdir = m.get_one::<String>("workdir").unwrap();
    let argv: Vec<&String> = m.get_many::<String>("vsearch").unwrap().collect();
    match rscompare::translate::translate(&argv, workdir) {
        Err(e) => {
            eprintln!("untranslatable: {e}");
            ExitCode::from(2)
        }
        Ok(p) if m.get_flag("json") => {
            println!("{}", serde_json::to_string_pretty(&p.to_json()).unwrap());
            ExitCode::SUCCESS
        }
        Ok(p) => {
            println!("{}", p.render());
            for o in &p.outputs {
                println!(
                    "# --{} {} <- {} ({})",
                    o.option, o.vsearch, o.rsearch, o.comparator
                );
            }
            for c in &p.caveats {
                println!("# caveat: {c}");
            }
            ExitCode::SUCCESS
        }
    }
}

fn list() -> ExitCode {
    for s in sources::all() {
        let seeds = (s.seeds)();
        let range = match (seeds.first(), seeds.last()) {
            (Some(a), Some(b)) if seeds.len() as u64 == b - a + 1 => format!("{a}-{b}"),
            (Some(_), _) => seeds
                .iter()
                .map(u64::to_string)
                .collect::<Vec<_>>()
                .join(","),
            _ => "-".into(),
        };
        println!("{:14} {range:10} {}", s.name, s.about);
    }
    ExitCode::SUCCESS
}

fn run(m: &clap::ArgMatches) -> ExitCode {
    let all: Vec<String> = sources::all().iter().map(|s| s.name.to_string()).collect();
    let mut names: Vec<String> = m.get_many::<String>("sources").unwrap().cloned().collect();
    if names == ["all"] {
        // every source but the 4 GB files of mkudb-long, and the datasets
        names = all
            .iter()
            .filter(|n| *n != "mkudb-long" && *n != "datasets")
            .cloned()
            .collect();
    }
    for n in &names {
        if !all.contains(n) {
            eprintln!("rscompare: unknown source {n}; see list");
            return ExitCode::from(2);
        }
    }
    let seeds = match m.get_one::<String>("seeds").map(|s| runner::parse_seeds(s)) {
        None => None,
        Some(Ok(s)) => Some(s),
        Some(Err(e)) => {
            eprintln!("rscompare: {e}");
            return ExitCode::from(2);
        }
    };
    let repo = runner::repo_root();
    let path = |name: &str, default: PathBuf| {
        let given = m.get_one::<String>(name).map_or(default, PathBuf::from);
        given.canonicalize().unwrap_or(given)
    };
    let opts = runner::Options {
        root: m
            .get_one::<String>("tmp")
            .map_or_else(runner::default_root, PathBuf::from),
        vsearch: path("vsearch", repo.join("compare/vsearch-build/vsearch")),
        bindir: path("bin", repo.join("target/release")),
        out: m.get_one::<String>("out").map(PathBuf::from),
        keep: m.get_flag("keep"),
        jobs: m
            .get_one::<usize>("jobs")
            .copied()
            .unwrap_or_else(runner::default_jobs),
    };
    match runner::run(&names, seeds.as_deref(), opts) {
        Err(e) => {
            eprintln!("rscompare: {e}");
            ExitCode::FAILURE
        }
        Ok(tally) => {
            runner::print_tally(&tally);
            let bad: u64 = tally
                .values()
                .map(|t| t.get("fail").copied().unwrap_or(0) + t.get("error").copied().unwrap_or(0))
                .sum();
            if bad > 0 {
                ExitCode::FAILURE
            } else {
                ExitCode::SUCCESS
            }
        }
    }
}

fn report(r: std::io::Result<()>) -> ExitCode {
    match r {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("rscompare: {e}");
            ExitCode::FAILURE
        }
    }
}
