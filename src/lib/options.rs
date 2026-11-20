//! Command-line options shared by the tools.

use clap::{Arg, ArgAction, ArgGroup, ArgMatches, Command, value_parser};

use crate::accept::Criteria;
use crate::scoring::Scoring;

use crate::label::{Label, Relabel};

/// adds the options that relabel sequences and add or remove label
/// annotations. Abundance and expected error annotations are left to each
/// tool, which decides whether it writes them.
pub fn label_args(command: Command) -> Command {
    let flag = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .action(ArgAction::SetTrue)
            .help(help)
    };
    command
        .arg(
            Arg::new("relabel")
                .long("relabel")
                .value_name("PREFIX")
                .help("Relabel sequences as PREFIX1, PREFIX2, ..."),
        )
        .arg(flag(
            "relabel-self",
            "Relabel sequences with the sequence itself",
        ))
        .arg(flag(
            "relabel-sha1",
            "Relabel sequences with the SHA1 digest of the sequence",
        ))
        .arg(flag(
            "relabel-md5",
            "Relabel sequences with the MD5 digest of the sequence",
        ))
        .group(ArgGroup::new("relabeling").args([
            "relabel",
            "relabel-self",
            "relabel-sha1",
            "relabel-md5",
        ]))
        .arg(flag(
            "relabel-keep",
            "Append the old label after the new one when relabelling",
        ))
        .arg(flag("xsize", "Remove abundance annotations (;size=)"))
        .arg(flag("xee", "Remove expected error annotations (;ee=)"))
        .arg(flag(
            "lengthout",
            "Annotate output labels with sequence length (;length=)",
        ))
        .arg(flag("xlength", "Remove length annotations (;length=)"))
        .arg(
            Arg::new("sample")
                .long("sample")
                .value_name("SAMPLE")
                .help("Annotate output labels with a sample name (;sample=)"),
        )
        .arg(
            Arg::new("label-suffix")
                .long("label-suffix")
                .value_name("SUFFIX")
                .help("Append this text to output labels"),
        )
}

/// returns the label settings of the options added by `label_args`, with
/// no abundance or expected error annotations.
pub fn label_from(matches: &ArgMatches) -> Label {
    let relabel = if let Some(prefix) = matches.get_one::<String>("relabel") {
        Relabel::Prefix(prefix.clone())
    } else if matches.get_flag("relabel-self") {
        Relabel::Sequence
    } else if matches.get_flag("relabel-sha1") {
        Relabel::Sha1
    } else if matches.get_flag("relabel-md5") {
        Relabel::Md5
    } else {
        Relabel::None
    };
    Label {
        relabel,
        keep: matches.get_flag("relabel-keep"),
        xsize: matches.get_flag("xsize"),
        xee: matches.get_flag("xee"),
        lengthout: matches.get_flag("lengthout"),
        xlength: matches.get_flag("xlength"),
        sample: matches.get_one::<String>("sample").cloned(),
        suffix: matches.get_one::<String>("label-suffix").cloned(),
        ..Label::default()
    }
}

/// a share from 0 to 1.
pub fn share(text: &str) -> Result<f64, String> {
    match text.parse::<f64>() {
        Ok(x) if (0.0..=1.0).contains(&x) => Ok(x),
        _ => Err("expected a number from 0 to 1".to_string()),
    }
}

/// the arguments of vsearch's acceptance criteria and scores,
/// shared by the tools that align.
pub fn search_args(command: Command) -> Command {
    let number = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("NUMBER")
            .value_parser(value_parser!(f64))
            .help(help)
    };
    let count = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .value_name("COUNT")
            .value_parser(value_parser!(i64).range(0..))
            .help(help)
    };
    let flag = |name: &'static str, help: &'static str| {
        Arg::new(name)
            .long(name)
            .action(ArgAction::SetTrue)
            .help(help)
    };
    let gaps = |text: &str, open: bool| -> Result<String, String> {
        Scoring::default()
            .set_gaps(text, open)
            .map(|()| text.to_string())
    };
    command
        .arg(
            Arg::new("id")
                .long("id")
                .value_name("SHARE")
                .value_parser(share)
                .help("Accept hits of at least this identity, 0 to 1"),
        )
        .arg(
            Arg::new("iddef")
                .long("iddef")
                .value_name("DEF")
                .value_parser(value_parser!(u8).range(0..=4))
                .default_value("2")
                .help("Identity definition, 0 to 4, as vsearch's"),
        )
        .arg(
            Arg::new("weak-id")
                .long("weak-id")
                .value_name("SHARE")
                .value_parser(share)
                .help("Weak hits: of at least this identity, below --id"),
        )
        .arg(count(
            "max-diffs",
            "Reject hits with more substitutions and interior gap columns",
        ))
        .arg(count("max-subs", "Reject hits with more substitutions"))
        .arg(count("max-gaps", "Reject hits with more interior gaps"))
        .arg(count("min-cols", "Reject hits with fewer interior columns"))
        .arg(flag(
            "left-just",
            "Reject hits with a terminal gap on the left",
        ))
        .arg(flag(
            "right-just",
            "Reject hits with a terminal gap on the right",
        ))
        .arg(number(
            "query-cov",
            "Reject hits covering less of the query, 0 to 1",
        ))
        .arg(number(
            "target-cov",
            "Reject hits covering less of the target, 0 to 1",
        ))
        .arg(number("max-id", "Reject hits of a higher identity, 0 to 1"))
        .arg(number(
            "mid",
            "Reject hits of a lower identity of letter pairs, in percent",
        ))
        .arg(count(
            "id-prefix",
            "Reject targets whose first residues differ from the query's",
        ))
        .arg(count(
            "id-suffix",
            "Reject targets whose last residues differ from the query's",
        ))
        .arg(number(
            "max-qt",
            "Reject targets unless the query's length is at most this times theirs",
        ))
        .arg(number(
            "min-qt",
            "Reject targets unless the query's length is at least this times theirs",
        ))
        .arg(number(
            "min-sl",
            "Reject pairs whose shorter length is below this share of the longer",
        ))
        .arg(number(
            "max-sl",
            "Reject pairs whose shorter length is above this share of the longer",
        ))
        .arg(count(
            "max-qsize",
            "Reject queries of an abundance above this",
        ))
        .arg(count(
            "min-tsize",
            "Reject targets of an abundance below this",
        ))
        .arg(number(
            "min-size-ratio",
            "Reject targets unless the query's abundance is at least this times theirs",
        ))
        .arg(number(
            "max-size-ratio",
            "Reject targets unless the query's abundance is at most this times theirs",
        ))
        .arg(flag("self", "Reject targets with the query's label"))
        .arg(flag("self-id", "Reject targets identical to the query"))
        .arg(
            Arg::new("match")
                .long("match")
                .value_name("SCORE")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i64).range(-32767..=32767))
                .default_value("2")
                .help("Score of a match"),
        )
        .arg(
            Arg::new("mismatch")
                .long("mismatch")
                .value_name("SCORE")
                .allow_negative_numbers(true)
                .value_parser(value_parser!(i64).range(-32767..=32767))
                .default_value("-4")
                .help("Score of a mismatch"),
        )
        .arg(
            Arg::new("gapopen")
                .long("gapopen")
                .value_name("PENALTIES")
                .value_parser(move |s: &str| gaps(s, true))
                .help("Gap opening penalties, as vsearch's, default to 20I/2E"),
        )
        .arg(
            Arg::new("gapext")
                .long("gapext")
                .value_name("PENALTIES")
                .value_parser(move |s: &str| gaps(s, false))
                .help("Gap extension penalties, as vsearch's, default to 2I/1E"),
        )
        .arg(flag("n-mismatch", "Count N as a mismatch"))
}

/// the arguments choosing the hits written of each query, for the tools
/// that write every hit.
pub fn hit_choice_args(command: Command) -> Command {
    command
        .arg(
            Arg::new("maxhits")
                .long("maxhits")
                .value_name("N")
                .value_parser(value_parser!(u64).range(1..))
                .help("Write at most this many hits of a query, default to all"),
        )
        .arg(
            Arg::new("top-hits-only")
                .long("top-hits-only")
                .action(ArgAction::SetTrue)
                .help("Write only the hits of the top identity"),
        )
}

/// the scores and criteria the arguments give.
pub fn search_options(matches: &ArgMatches) -> (Scoring, Criteria) {
    let mut scoring = Scoring {
        match_score: *matches.get_one::<i64>("match").unwrap(),
        mismatch: *matches.get_one::<i64>("mismatch").unwrap(),
        n_mismatch: matches.get_flag("n-mismatch"),
        ..Scoring::default()
    };
    if let Some(text) = matches.get_one::<String>("gapopen") {
        scoring.set_gaps(text, true).unwrap();
    }
    if let Some(text) = matches.get_one::<String>("gapext") {
        scoring.set_gaps(text, false).unwrap();
    }
    let mut c = Criteria::default();
    let f = |name: &str| matches.get_one::<f64>(name).copied();
    let i = |name: &str| matches.get_one::<i64>(name).copied();
    if let Some(v) = f("id") {
        c.id = v;
    }
    if let Some(v) = f("weak-id") {
        c.weak_id = v;
    }
    // as vsearch's, weak hits are never above --id
    if c.weak_id > c.id {
        c.weak_id = c.id;
    }
    c.iddef = usize::from(*matches.get_one::<u8>("iddef").unwrap());
    c.max_diffs = i("max-diffs").unwrap_or(c.max_diffs);
    c.max_subs = i("max-subs").unwrap_or(c.max_subs);
    c.max_gaps = i("max-gaps").unwrap_or(c.max_gaps);
    c.min_cols = i("min-cols").unwrap_or(c.min_cols);
    c.left_just = matches.get_flag("left-just");
    c.right_just = matches.get_flag("right-just");
    c.query_cov = f("query-cov").unwrap_or(c.query_cov);
    c.target_cov = f("target-cov").unwrap_or(c.target_cov);
    c.max_id = f("max-id").unwrap_or(c.max_id);
    c.mid = f("mid").unwrap_or(c.mid);
    c.id_prefix = i("id-prefix").map_or(0, |v| v as usize);
    c.id_suffix = i("id-suffix").map_or(0, |v| v as usize);
    c.max_qt = f("max-qt").unwrap_or(c.max_qt);
    c.min_qt = f("min-qt").unwrap_or(c.min_qt);
    c.min_sl = f("min-sl").unwrap_or(c.min_sl);
    c.max_sl = f("max-sl").unwrap_or(c.max_sl);
    c.sizes.max_query_size = i("max-qsize").unwrap_or(c.sizes.max_query_size);
    c.sizes.min_target_size = i("min-tsize").unwrap_or(c.sizes.min_target_size);
    c.sizes.min_size_ratio = f("min-size-ratio").unwrap_or(c.sizes.min_size_ratio);
    c.sizes.max_size_ratio = f("max-size-ratio").unwrap_or(c.sizes.max_size_ratio);
    c.self_label = matches.get_flag("self");
    c.self_id = matches.get_flag("self-id");
    (scoring, c)
}
