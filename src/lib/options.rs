//! Command-line options shared by the tools.

use clap::{Arg, ArgAction, ArgGroup, ArgMatches, Command};

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
