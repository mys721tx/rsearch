//! vsearch's options and a parser for its command lines.
//!
//! The option names, and whether each takes a value, come from the option
//! table in vsearch's `src/cli.cc` (`../vsearch`, e1551256).

use std::collections::{HashMap, HashSet};
use std::fmt;

/// Options that take no value.
pub const FLAGS: &[&str] = &[
    "acceptall",
    "allow_fewer",
    "bzip2_decompress",
    "centroid_sizeout",
    "clusterout_id",
    "clusterout_sort",
    "cons_truncate",
    "eeout",
    "fasta_score",
    "fastq_allowmergestagger",
    "fastq_eeout",
    "fastq_nostagger",
    "fastq_qout_max",
    "fastq_solexa",
    "fulldp",
    "gzip_decompress",
    "h",
    "hardmask",
    "help",
    "label_substr_match",
    "leftjust",
    "lengthout",
    "n_mismatch",
    "no_progress",
    "notrunclabels",
    "output_no_hits",
    "quiet",
    "relabel_keep",
    "relabel_md5",
    "relabel_self",
    "relabel_sha1",
    "rightjust",
    "samheader",
    "self",
    "selfid",
    "sff_clip",
    "sintax_random",
    "sizein",
    "sizeorder",
    "sizeout",
    "top_hits_only",
    "uc_allhits",
    "uchimeout5",
    "usersort",
    "v",
    "version",
    "xee",
    "xlength",
    "xsize",
];

/// Options that take a value, the commands included.
pub const VALUED: &[&str] = &[
    "abskew",
    "alignwidth",
    "allpairs_global",
    "alnout",
    "band",
    "biomout",
    "blast6out",
    "borderline",
    "centroids",
    "chimeras",
    "chimeras_denovo",
    "chimeras_diff_pct",
    "chimeras_length_min",
    "chimeras_parents_max",
    "chimeras_parts",
    "cluster_fast",
    "cluster_size",
    "cluster_smallmem",
    "cluster_unoise",
    "clusters",
    "consout",
    "cut",
    "cut_pattern",
    "db",
    "dbmask",
    "dbmatched",
    "dbnotmatched",
    "derep_fulllength",
    "derep_id",
    "derep_prefix",
    "derep_smallmem",
    "dn",
    "ee_cutoffs",
    "eetabbedout",
    "fasta2fastq",
    "fasta_width",
    "fastaout",
    "fastaout_discarded",
    "fastaout_discarded_rev",
    "fastaout_notmerged_fwd",
    "fastaout_notmerged_rev",
    "fastaout_orphans",
    "fastaout_orphans_rev",
    "fastaout_rev",
    "fastapairs",
    "fastq_ascii",
    "fastq_asciiout",
    "fastq_chars",
    "fastq_convert",
    "fastq_eestats",
    "fastq_eestats2",
    "fastq_filter",
    "fastq_join",
    "fastq_maxdiffpct",
    "fastq_maxdiffs",
    "fastq_maxee",
    "fastq_maxee_rate",
    "fastq_maxlen",
    "fastq_maxmergelen",
    "fastq_maxns",
    "fastq_mergepairs",
    "fastq_minlen",
    "fastq_minmergelen",
    "fastq_minovlen",
    "fastq_minqual",
    "fastq_qmax",
    "fastq_qmaxout",
    "fastq_qmin",
    "fastq_qminout",
    "fastq_stats",
    "fastq_stripleft",
    "fastq_stripright",
    "fastq_tail",
    "fastq_truncee",
    "fastq_truncee_rate",
    "fastq_trunclen",
    "fastq_trunclen_keep",
    "fastq_truncqual",
    "fastqout",
    "fastqout_discarded",
    "fastqout_discarded_rev",
    "fastqout_notmerged_fwd",
    "fastqout_notmerged_rev",
    "fastqout_orphans",
    "fastqout_orphans_rev",
    "fastqout_rev",
    "fastx_filter",
    "fastx_getseq",
    "fastx_getseqs",
    "fastx_getsubseq",
    "fastx_mask",
    "fastx_revcomp",
    "fastx_subsample",
    "fastx_syncpairs",
    "fastx_uniques",
    "gapext",
    "gapopen",
    "hspw",
    "id",
    "iddef",
    "idprefix",
    "idsuffix",
    "join_padgap",
    "join_padgapq",
    "label",
    "label_field",
    "label_suffix",
    "label_word",
    "label_words",
    "labels",
    "lca_cutoff",
    "lcaout",
    "length_cutoffs",
    "log",
    "makeudb_usearch",
    "maskfasta",
    "match",
    "matched",
    "max_unmasked_pct",
    "maxaccepts",
    "maxdiffs",
    "maxgaps",
    "maxhits",
    "maxid",
    "maxqsize",
    "maxqt",
    "maxrejects",
    "maxseqlength",
    "maxsize",
    "maxsizeratio",
    "maxsl",
    "maxsubs",
    "maxuniquesize",
    "mid",
    "min_unmasked_pct",
    "mincols",
    "mindiffs",
    "mindiv",
    "minh",
    "minhsp",
    "minqt",
    "minseqlength",
    "minsize",
    "minsizeratio",
    "minsl",
    "mintsize",
    "minuniquesize",
    "minwordmatches",
    "mismatch",
    "mothur_shared_out",
    "msaout",
    "nonchimeras",
    "notmatched",
    "notmatchedfq",
    "orient",
    "otutabout",
    "output",
    "pattern",
    "profile",
    "qmask",
    "qsegout",
    "query_cov",
    "randseed",
    "read_separators",
    "relabel",
    "rereplicate",
    "reverse",
    "rowlen",
    "samout",
    "sample",
    "sample_pct",
    "sample_size",
    "scramble",
    "scramble_kmer",
    "search_exact",
    "sff_convert",
    "shuffle",
    "sintax",
    "sintax_cutoff",
    "slots",
    "sortbylength",
    "sortbysize",
    "strand",
    "subseq_end",
    "subseq_start",
    "tabbedout",
    "target_cov",
    "threads",
    "topn",
    "tsegout",
    "uc",
    "uchime2_denovo",
    "uchime3_denovo",
    "uchime_denovo",
    "uchime_ref",
    "uchimealns",
    "uchimeout",
    "udb2fasta",
    "udbinfo",
    "udbstats",
    "unoise_alpha",
    "usearch_global",
    "userfields",
    "userout",
    "weak_id",
    "wordlength",
    "xdrop_nw",
    "xn",
];

/// The 48 commands, and `--help`, `--version`, `-h`, `-v`.
pub const COMMANDS: &[&str] = &[
    "allpairs_global",
    "chimeras_denovo",
    "cluster_fast",
    "cluster_size",
    "cluster_smallmem",
    "cluster_unoise",
    "cut",
    "derep_fulllength",
    "derep_id",
    "derep_prefix",
    "derep_smallmem",
    "fasta2fastq",
    "fastq_chars",
    "fastq_convert",
    "fastq_eestats",
    "fastq_eestats2",
    "fastq_filter",
    "fastq_join",
    "fastq_mergepairs",
    "fastq_stats",
    "fastx_filter",
    "fastx_getseq",
    "fastx_getseqs",
    "fastx_getsubseq",
    "fastx_mask",
    "fastx_revcomp",
    "fastx_subsample",
    "fastx_syncpairs",
    "fastx_uniques",
    "makeudb_usearch",
    "maskfasta",
    "orient",
    "rereplicate",
    "scramble",
    "search_exact",
    "sff_convert",
    "shuffle",
    "sintax",
    "sortbylength",
    "sortbysize",
    "uchime2_denovo",
    "uchime3_denovo",
    "uchime_denovo",
    "uchime_ref",
    "udb2fasta",
    "udbinfo",
    "udbstats",
    "usearch_global",
    "help",
    "version",
    "h",
    "v",
];

const INFO: &[&str] = &["help", "version", "h", "v"];

/// A vsearch command line the translator cannot turn into rsearch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranslationError(pub String);

impl fmt::Display for TranslationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for TranslationError {}

/// `Err(TranslationError)` from a message.
pub fn err<T>(message: impl Into<String>) -> Result<T, TranslationError> {
    Err(TranslationError(message.into()))
}

/// A value of an option: a flag, or the text given.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Value {
    Flag,
    Text(String),
}

/// A parsed vsearch command line: the command, its argument, and the other
/// options in order. The translator takes options as it uses them; what is
/// left untaken is untranslatable.
#[derive(Debug, Clone)]
pub struct Args {
    pub command: String,
    pub input: Option<String>,
    pub opts: HashMap<String, Value>,
    pub order: Vec<String>,
    pub used: HashSet<String>,
}

impl Args {
    /// Parses vsearch's argv, with or without the program name.
    pub fn parse<S: AsRef<str>>(argv: &[S]) -> Result<Args, TranslationError> {
        let mut words: Vec<&str> = argv.iter().map(AsRef::as_ref).collect();
        if words.first().is_some_and(|w| !w.starts_with('-')) {
            words.remove(0);
        }
        let mut args = Args {
            command: String::new(),
            input: None,
            opts: HashMap::new(),
            order: Vec::new(),
            used: HashSet::new(),
        };
        let mut i = 0;
        while i < words.len() {
            let tok = words[i];
            if !tok.starts_with('-') || tok == "-" || tok == "--" {
                return err(format!("unexpected argument {}", quote(tok)));
            }
            let mut name = tok.trim_start_matches('-').to_string();
            let mut value = None;
            if let Some(k) = name.find('=') {
                value = Some(name[k + 1..].to_string());
                name.truncate(k);
            }
            let is_command = COMMANDS.contains(&name.as_str());
            let value = if FLAGS.contains(&name.as_str()) && !is_command {
                if value.is_some() {
                    return err(format!("--{name} takes no value"));
                }
                Value::Flag
            } else if VALUED.contains(&name.as_str()) || is_command {
                match value {
                    Some(v) => Value::Text(v),
                    None if INFO.contains(&name.as_str()) => Value::Flag,
                    None => {
                        i += 1;
                        match words.get(i) {
                            Some(v) => Value::Text(v.to_string()),
                            None => return err(format!("--{name} needs a value")),
                        }
                    }
                }
            } else {
                return err(format!("unknown vsearch option --{name}"));
            };
            if is_command {
                if !args.command.is_empty() {
                    return err(format!("two commands: --{} and --{name}", args.command));
                }
                args.command = name;
                if let Value::Text(v) = value {
                    args.input = Some(v);
                }
            } else {
                if !args.opts.contains_key(&name) {
                    args.order.push(name.clone());
                }
                args.opts.insert(name, value);
            }
            i += 1;
        }
        if args.command.is_empty() {
            return err("no vsearch command");
        }
        Ok(args)
    }

    pub fn has(&self, name: &str) -> bool {
        self.opts.contains_key(name)
    }

    /// The option's text (marking it used); a flag gives "".
    pub fn take(&mut self, name: &str) -> Option<String> {
        let v = self.opts.get(name)?.clone();
        self.used.insert(name.to_string());
        Some(match v {
            Value::Flag => String::new(),
            Value::Text(t) => t,
        })
    }

    /// The option's text without marking it used.
    pub fn peek(&self, name: &str) -> Option<&str> {
        match self.opts.get(name)? {
            Value::Flag => Some(""),
            Value::Text(t) => Some(t),
        }
    }

    /// Whether a flag is given (marking it used).
    pub fn flag(&mut self, name: &str) -> bool {
        self.take(name).is_some()
    }

    /// The option's value as a number, or the default.
    pub fn take_num(&mut self, name: &str, default: f64) -> Result<f64, TranslationError> {
        match self.take(name) {
            None => Ok(default),
            Some(v) => v
                .parse::<f64>()
                .or_else(|_| err(format!("--{name} {v}: not a number"))),
        }
    }

    /// Replaces an option's value (a decompressed input, say).
    pub fn set(&mut self, name: &str, value: String) {
        self.opts.insert(name.to_string(), Value::Text(value));
    }

    /// The options never taken, in the order given.
    pub fn untaken(&self) -> Vec<String> {
        self.order
            .iter()
            .filter(|n| !self.used.contains(*n))
            .cloned()
            .collect()
    }
}

/// Python's repr of a string, which the error messages use.
pub fn quote(s: &str) -> String {
    if s.contains('\'') && !s.contains('"') {
        format!("\"{s}\"")
    } else {
        format!("'{}'", s.replace('\'', "\\'"))
    }
}
