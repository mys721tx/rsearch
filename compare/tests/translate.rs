//! The translator on command lines: the pipelines of the vsearch commands,
//! the option names, and the errors and caveats.
//!
//! Files named here do not exist. No input is therefore known to be FASTA,
//! and FASTA outputs always go through conv.

use std::collections::HashMap;
use std::path::Path;

use rscompare::options::{Args, COMMANDS, Value};
use rscompare::translate::{Pipeline, THREADED_TOOLS, translate};

const W: &str = "w";

fn words(cmd: &str) -> Vec<&str> {
    cmd.split_whitespace().collect()
}

fn pipe(cmd: &str) -> Pipeline {
    translate(&words(cmd), W).unwrap_or_else(|e| panic!("{cmd}: {e}"))
}

fn render(cmd: &str) -> String {
    pipe(cmd).render()
}

fn outputs(cmd: &str) -> HashMap<String, (String, String)> {
    pipe(cmd)
        .outputs
        .into_iter()
        .map(|o| (o.option, (o.rsearch, o.comparator)))
        .collect()
}

fn caveats(cmd: &str) -> Vec<String> {
    pipe(cmd).caveats
}

#[track_caller]
fn has(haystack: &str, needle: &str) {
    assert!(haystack.contains(needle), "{needle:?} not in\n{haystack}");
}

#[track_caller]
fn bad(cmd: &str, text: &str) {
    match translate(&words(cmd), W) {
        Ok(p) => panic!("{cmd} translated:\n{}", p.render()),
        Err(e) => assert!(e.0.contains(text), "{text:?} not in {:?}", e.0),
    }
}

// the preprocessing pipelines, one by one

#[test]
fn filter_trunclen() {
    assert_eq!(
        render("--fastx_filter in --fastq_trunclen 200 --fastq_maxee 1 --relabel R --fastaout out"),
        "sqtrim -i in --len 200 | sqfilt --min-len 200 --max-ee 1 | relab --relabel R | conv \
         --fasta > w/tmp4.fastaout"
    );
}

#[test]
fn filter_pairs() {
    assert_eq!(
        render(
            "--fastq_filter in --reverse rev --fastq_stripleft 3 --fastq_maxee 1 --fastqout k \
             --fastqout_rev k2"
        ),
        "sqtrim -i in --strip-left 3 > w/tmp1.sqtrim\nsqtrim -i rev --strip-left 3 > \
         w/tmp2.sqtrim\nsqfilt --reverse w/tmp2.sqtrim --rev-out w/tmp4.rev_kept --max-ee 1 < \
         w/tmp1.sqtrim > w/tmp3.kept"
    );
}

#[test]
fn uniques() {
    assert_eq!(
        render(
            "--fastx_uniques in --minseqlength 32 --minuniquesize 2 --topn 10 --relabel U \
             --sizeout --fastaout o"
        ),
        "sqfilt -i in --min-len 32 --max-len 50000 | derep --no-sizein | sqfilt --min-size 2 | \
         topn 10 | relab --relabel U --sizeout | conv --fasta > w/tmp6.fastaout"
    );
}

#[test]
fn derep_smallmem() {
    assert_eq!(
        render("--derep_smallmem in --fastaout out"),
        "sqfilt -i in --max-len 50000 | derep --unsorted --no-sizein --no-sizeout | conv \
         --fasta > w/tmp3.fastaout"
    );
}

#[test]
fn rereplicate() {
    assert_eq!(
        render("--rereplicate in --sizeout --output out"),
        "rerep -i in | conv --fasta > w/tmp2.output"
    );
    assert_eq!(
        render("--rereplicate in --output out"),
        "rerep -i in --no-sizeout | conv --fasta > w/tmp2.output"
    );
    // relab drops ;size= unless it writes it again
    assert_eq!(
        render("--rereplicate in --relabel R --sizeout --output out"),
        "rerep -i in | relab --relabel R --sizeout | conv --fasta > w/tmp3.output"
    );
    assert_eq!(
        render("--rereplicate in --relabel R --output out"),
        "rerep -i in --no-sizeout | relab --relabel R | conv --fasta > w/tmp3.output"
    );
}

#[test]
fn join() {
    assert_eq!(
        render("--fastq_join in --reverse rev --relabel J --fastaout out"),
        "joinpa -i in --reverse rev | relab --relabel J | conv --fasta > w/tmp3.fastaout"
    );
}

#[test]
fn mergepairs_filters() {
    assert_eq!(
        render(
            "--fastq_mergepairs in --reverse rev --fastq_truncqual 2 --fastq_minlen 50 \
             --fastq_maxee 1 --fastqout m"
        ),
        "sqtrim -i in --qual 2 > w/tmp1.sqtrim\nsqtrim -i rev --qual 2 > w/tmp2.sqtrim\nsqfilt \
         --reverse w/tmp2.sqtrim --rev-out w/tmp3.rev --min-len 50 < w/tmp1.sqtrim | mergpa \
         --reverse w/tmp3.rev | sqfilt --max-ee 1 > w/tmp6.sqfilt"
    );
}

#[test]
fn mergepairs_lengths() {
    assert_eq!(
        render(
            "--fastq_mergepairs in --reverse rev --fastq_minmergelen 100 --fastq_maxmergelen 300 \
             --fastqout m"
        ),
        "mergpa -i in --reverse rev | sqfilt --min-len 100 --max-len 300 > w/tmp2.sqfilt"
    );
}

#[test]
fn getsubseq() {
    assert_eq!(
        render("--fastx_getsubseq in --label s1 --subseq_start 3 --subseq_end 6 --fastaout o"),
        "getseq -i in --label s1 | subseq --start 3 --end 6 | conv --fasta > w/tmp3.fastaout"
    );
}

#[test]
fn sortbysize() {
    assert_eq!(
        render("--sortbysize in --minsize 2 --topn 5 --output out"),
        "sqfilt -i in --max-len 50000 --min-size 2 | sortby --size | topn 5 | conv --fasta > \
         w/tmp4.output"
    );
}

#[test]
fn fastq_convert() {
    assert_eq!(
        render("--fastq_convert in --fastq_ascii 64 --relabel R --fastqout o"),
        "conv -i in --ascii 64 | relab --relabel R > w/tmp2.relab"
    );
}

#[test]
fn cut() {
    assert_eq!(
        render("--cut in --cut_pattern G^AATT_C --fastaout f --fastaout_discarded d"),
        "redig -i in --pattern 'G^AATT_C' --uncut w/tmp2.uncut | conv --fasta > \
         w/tmp3.fastaout\nconv --fasta < w/tmp2.uncut > w/tmp4.fastaout_discarded"
    );
}

#[test]
fn mask() {
    assert_eq!(
        render("--fastx_mask in --fastaout out"),
        "lcmask -i in | sqfilt --max-len 50000 | conv --fasta > w/tmp3.fastaout"
    );
    assert_eq!(
        render("--maskfasta in --output out"),
        "lcmask -i in | sqfilt --max-len 50000 | conv --fasta > w/tmp3.output"
    );
    has(
        &render("--fastx_mask in --qmask soft --hardmask --fastaout o"),
        "lcmask -i in --hard --no-dust |",
    );
    has(
        &render("--fastx_mask in --min_unmasked_pct 50 --fastaout o"),
        "lcmask -i in | sqfilt --max-len 50000 --min-upper-pct 50",
    );
    has(
        &render("--fastx_mask in --hardmask --min_unmasked_pct 50 --fastaout o"),
        "lcmask -i in --hard | sqfilt --max-len 50000 --min-non-n-pct 50",
    );
}

// the analysis pipelines, one by one

#[test]
fn cluster_fast() {
    assert_eq!(
        render("--cluster_fast in --id 0.97 --sizein --sizeout --centroids c --uc u"),
        "sqfilt -i in --min-len 32 --max-len 50000 | lcmask | sortby --length | clust \
         --maxrejects 8 --id 0.97 --hits w/tmp4.uc | conv --fasta > w/tmp6.centroids"
    );
    assert_eq!(
        outputs("--cluster_fast in --id 0.97 --uc u")["uc"].1,
        "uc-fast"
    );
}

#[test]
fn cluster_size() {
    has(
        &render("--cluster_size in --id 0.97 --sizein --sizeout --centroids c"),
        "sortby --size | clust --id 0.97 | conv --fasta",
    );
}

#[test]
fn cluster_smallmem() {
    has(
        &render("--cluster_smallmem in --usersort --id 0.97 --centroids c"),
        "lcmask | clust --id 0.97 --no-sizein --no-sizeout",
    );
}

#[test]
fn cluster_unoise() {
    has(
        &render("--cluster_unoise in --sizein --sizeout --centroids c"),
        "sqfilt -i in --min-len 32 --max-len 50000 --min-size 8 | lcmask | sortby --size | \
         unoise | conv --fasta",
    );
    let r = render(
        "--cluster_unoise in --minsize 2 --weak_id 0.8 --unoise_alpha 3 --id 0.5 --centroids c",
    );
    has(&r, "unoise --id 0.8 --alpha 3 --no-sizein");
    assert!(!r.contains("0.5"));
}

#[test]
fn cluster_outputs() {
    assert_eq!(
        render(
            "--cluster_fast in --id 0.97 --sizein --sizeout --consout cons --clusters cl \
             --fasta_width 0"
        ),
        "sqfilt -i in --min-len 32 --max-len 50000 | lcmask | sortby --length > \
         w/tmp3.sortby\nclust --maxrejects 8 --id 0.97 --hits w/tmp4.uc < w/tmp3.sortby > \
         w/tmp5.centroids\nclspl --hits w/tmp4.uc -i w/tmp3.sortby --prefix w/clusters > \
         w/tmp6.clspl\nclmsa --hits w/tmp4.uc -i w/tmp3.sortby --consensus > w/tmp7.consout"
    );
}

#[test]
fn usearch_global_alnout() {
    assert_eq!(
        render("--usearch_global q --db db --id 0.97 --alnout a"),
        "lcmask -i q > w/tmp1.lcmask\nlcmask -i db | sqfilt --min-len 32 --max-len 50000 > \
         w/tmp3.db\ndbsrch --db w/tmp3.db --id 0.97 < w/tmp1.lcmask > w/tmp4.uc\nhitfmt \
         w/tmp4.uc -i w/tmp1.lcmask --db w/tmp3.db --alignments > w/tmp5.alnout"
    );
}

#[test]
fn usearch_global_tables() {
    let r = render("--usearch_global q --db db --id 0.97 --otutabout t --lcaout l");
    has(&r, "otutab w/tmp4.uc --db w/tmp3.db > w/tmp5.otutabout");
    has(&r, "lcatax w/tmp4.uc --db w/tmp3.db > w/tmp6.lcaout");
    assert_eq!(
        outputs("--usearch_global q --db db --id 1 --biomout b")["biomout"].1,
        "biom"
    );
}

#[test]
fn usearch_global_matched() {
    let r = render("--usearch_global q --db db --id 0.97 --matched m --dbnotmatched n");
    has(
        &r,
        "getseq --hits w/tmp4.uc -i w/tmp1.lcmask | conv --fasta",
    );
    has(
        &r,
        "getseq --hits w/tmp4.uc -i w/tmp3.db --targets --invert",
    );
}

#[test]
fn usearch_global_qmask_none() {
    let r = render(
        "--usearch_global q --db db --id 0.97 --qmask none --alnout a --rowlen 30 --gapopen \
         20I/2E --n_mismatch",
    );
    has(&r, "lcmask -i q --clear");
    // hitfmt shows the queries as they are, with the search's scores
    has(
        &r,
        "hitfmt w/tmp4.uc -i q --db w/tmp3.db --alignments --rowlen 30 --gapopen 20I/2E \
         --n-mismatch",
    );
}

#[test]
fn uc_top_hit_and_dbmatched() {
    let r = render("--usearch_global q --db db --id 0.97 --uc u");
    has(&r, "--maxhits 1 < w/tmp1.lcmask > w/tmp4.uc");
    let r = render(
        "--usearch_global q --db db --id 0.97 --dbmatched m --maxhits 2 --uc u --uc_allhits",
    );
    // --dbmatched counts every accepted hit: getseq --targets reads a search
    // without --maxhits
    has(&r, "--id 0.97 --maxhits 2 < w/tmp1.lcmask > w/tmp4.uc");
    has(&r, "--id 0.97 < w/tmp1.lcmask > w/tmp5.uc");
    has(&r, "getseq --hits w/tmp5.uc -i w/tmp3.db --targets");
}

#[test]
fn orient() {
    assert_eq!(
        render("--orient q --db db --fastaout o"),
        "sqfilt -i db --max-len 50000 | lcmask > w/tmp2.lcmask\ndbornt -i q --db \
         w/tmp2.lcmask | conv --fasta > w/tmp4.fastaout"
    );
    assert!(!render("--orient q --db db --dbmask soft --fastaout o").contains("lcmask"));
}

#[test]
fn uchime() {
    assert_eq!(
        render("--uchime_denovo in --nonchimeras n"),
        "sqfilt -i in --max-len 50000 | lcmask | sortby --size | uchime | conv --fasta > \
         w/tmp5.nonchimeras"
    );
    has(
        &render("--uchime3_denovo in --nonchimeras n --uchimeout t"),
        "uchime --uchime2 --abskew 16 --table w/uchimeout",
    );
    assert_eq!(
        render("--uchime_ref in --db ref --nonchimeras n"),
        "lcmask -i ref | sqfilt --max-len 50000 > w/tmp2.db\nuchime -i in --db w/tmp2.db | \
         conv --fasta > w/tmp4.nonchimeras"
    );
    // vsearch masks no query of --uchime_ref: dust and soft both skip the
    // words of its lower case as read, and the outputs are the queries as
    // read
    assert_eq!(
        render("--uchime_ref in --db ref --qmask soft --dbmask soft --nonchimeras n"),
        "sqfilt -i ref --max-len 50000 > w/tmp1.db\nuchime -i in --db w/tmp1.db | conv --fasta \
         > w/tmp3.nonchimeras"
    );
    assert!(
        caveats("--uchime_ref in --db ref --qmask none --nonchimeras n --threads 1")
            .iter()
            .any(|c| c.contains("--qmask none"))
    );
}

#[test]
fn makeudb() {
    assert_eq!(
        render("--makeudb_usearch in --output db.udb"),
        "sqfilt -i in --min-len 32 --max-len 50000 | lcmask | mkudb > w/tmp3.udb"
    );
    has(
        &render("--makeudb_usearch in --dbmask none --output o"),
        "| lcmask --clear | mkudb",
    );
    has(
        &render("--makeudb_usearch in --dbmask soft --output o"),
        "| lcmask --clear | mkudb",
    );
    has(
        &render("--makeudb_usearch in --hardmask --output o"),
        "| lcmask --clear | lcmask --hard | mkudb",
    );
}

#[test]
fn sintax() {
    assert_eq!(
        render("--sintax q --db db --randseed 1 --tabbedout t"),
        "sqfilt -i db --min-len 32 --max-len 50000 > w/tmp1.db\nlcmask -i q --clear | sintax \
         --db w/tmp1.db --seed=1 > w/tmp3.sintax"
    );
}

#[test]
fn search_exact() {
    assert_eq!(
        render("--search_exact q --db db --uc u --uc_allhits"),
        "sqfilt -i db --max-len 50000 > w/tmp1.db\nexsrch -i q --db w/tmp1.db > w/tmp2.uc"
    );
    has(
        &render("--search_exact q --db db --uc u"),
        "exsrch -i q --db w/tmp1.db --maxhits 1",
    );
    // hitfmt scores exact hits as vsearch's exact search does
    has(
        &render(
            "--search_exact q --db db --qmask none --dbmask none --match 3 --userout u --userfields raw",
        ),
        "--userout raw --match 3 --exact",
    );
}

#[test]
fn search_hits_tie() {
    // hits tying on vsearch's sort key may come in either order. Their
    // comparators read the queries
    for cmd in [
        "--search_exact q --db db --strand both",
        "--usearch_global q --db db --id 0.9 --strand both",
    ] {
        let o = outputs(&format!("{cmd} --uc u --blast6out b --samout s"));
        assert_eq!(o["uc"].1, "uc-ties:q", "{cmd}");
        assert_eq!(o["blast6out"].1, "blast6-ties:q", "{cmd}");
        assert_eq!(o["samout"].1, "sam-ties:q", "{cmd}");
    }
    // the other formats stay exact
    let o = outputs("--search_exact q --db db --userout u --userfields query+target");
    assert_eq!(o["userout"].1, "bytes");
}

#[test]
fn allpairs() {
    assert_eq!(
        render("--allpairs_global in --id 0.9 --uc u"),
        "sqfilt -i in --max-len 50000 | lcmask | allpr --id 0.9 > w/tmp3.uc"
    );
}

#[test]
fn udb2fasta() {
    assert_eq!(
        render("--udb2fasta db.udb --output out"),
        "conv -i db.udb --fasta > w/tmp1.output"
    );
}

// options

#[test]
fn names() {
    let r = render(
        "--allpairs_global in --qmask soft --minseqlength 1 --maxseqlength 100 --id 0.5 \
         --weak_id 0.4 --maxdiffs 3 --maxsubs 2 --maxgaps 1 --mincols 4 --query_cov 0.5 \
         --target_cov 0.6 --maxid 0.99 --mid 80 --idprefix 2 --idsuffix 3 --maxqt 2 --minqt \
         0.5 --minsl 0.5 --maxsl 2 --maxqsize 9 --mintsize 2 --minsizeratio 0.5 --maxsizeratio \
         3 --leftjust --rightjust --self --selfid --n_mismatch --iddef 1 --match 3 --mismatch -2 \
         --gapopen 20 --gapext 2 --maxhits 4 --top_hits_only --maxaccepts 3 --wordlength 9 --uc u",
    );
    assert_eq!(
        r,
        "sqfilt -i in --max-len 100 | allpr --id 0.5 --iddef 1 --weak-id 0.4 --max-diffs 3 \
         --max-subs 2 --max-gaps 1 --min-cols 4 --query-cov 0.5 --target-cov 0.6 --max-id 0.99 \
         --mid 80 --id-prefix 2 --id-suffix 3 --max-qt 2 --min-qt 0.5 --min-sl 0.5 --max-sl 2 \
         --max-qsize 9 --min-tsize 2 --min-size-ratio 0.5 --max-size-ratio 3 --match 3 \
         --mismatch -2 --gapopen 20 --gapext 2 --left-just --right-just --self --self-id \
         --n-mismatch --maxhits 4 --top-hits-only > w/tmp2.uc"
    );
}

#[test]
fn word_options() {
    has(
        &render(
            "--usearch_global q --db db --qmask soft --dbmask soft --id 1 --strand both \
             --maxaccepts 0 --maxrejects 0 --wordlength 10 --minwordmatches 3 --uc u",
        ),
        "dbsrch -i q --db w/tmp1.db --id 1 --strand both --maxaccepts 0 --maxrejects 0 \
         --word-len 10 --min-word-matches 3",
    );
}

#[test]
fn quality() {
    assert_eq!(
        render(
            "--fastq_convert in --fastq_ascii 64 --fastq_qmin 2 --fastq_qmax 40 --fastq_asciiout \
             33 --fastq_qminout 3 --fastq_qmaxout 39 --fastq_solexa --fastqout o"
        ),
        "conv -i in --ascii 64 --qmin 2 --qmax 40 --ascii-out 33 --qmin-out 3 --qmax-out 39 \
         --solexa > w/tmp1.conv"
    );
}

#[test]
fn labels() {
    has(
        &render(
            "--fastq_convert in --fastq_ascii 64 --relabel R --sample S --label_suffix ;x \
             --relabel_keep --sizeout --xsize --lengthout --xlength --xee --eeout --fastqout o",
        ),
        "relab --relabel R --sample S --label-suffix ';x' --relabel-keep --sizeout --xsize \
         --lengthout --xlength --xee --ee --ascii 64",
    );
}

#[test]
fn fasta_width() {
    has(
        &render("--fastx_revcomp in --fasta_width 0 --fastaout o"),
        "conv --fasta --fasta-width 0",
    );
    has(
        &render("--fastx_revcomp in --fasta_width 7 --fastaout o"),
        "conv --fasta --fasta-width 7",
    );
}

#[test]
fn seed() {
    has(
        &render("--shuffle in --randseed -3 --output o"),
        "--seed=-3",
    );
    has(
        &caveats("--shuffle in --output o").join(" "),
        "no --randseed",
    );
}

#[test]
fn subsample() {
    has(
        &render(
            "--fastx_subsample in --sample_size 5 --sizein --sizeout --randseed 1 --fastqout k \
             --fastqout_discarded d",
        ),
        "subsmp --count 5 --by-size --seed=1 --discarded",
    );
}

#[test]
fn eestats2() {
    assert_eq!(
        render("--fastq_eestats2 in --length_cutoffs 3,*,7 --ee_cutoffs 0.5,1 --output o"),
        "eegrid -i in --min-len 3 --len-step 7 --max-ee 0.5,1 > w/tmp1.eegrid"
    );
}

#[test]
fn logs() {
    let o = |cmd: &str, opt: &str| outputs(cmd)[opt].clone();
    assert_eq!(
        outputs("--fastq_chars in --fastq_tail 3 --log l"),
        HashMap::from([(
            "log".to_string(),
            ("w/log".to_string(), "log-report".to_string())
        )])
    );
    assert_eq!(o("--fastq_stats in --log l", "log").1, "log-report");
    assert_eq!(o("--udbinfo d --log l", "log").1, "udbinfo-log");
    assert_eq!(o("--udbstats d --log l", "log").1, "udbstats-log");
    assert_eq!(
        o("--orient q --db d --tabbedout t --log l", "log").1,
        "orient-log"
    );
    assert_eq!(
        o("--sintax q --db d --tabbedout t --log l", "log").1,
        "sintax-log"
    );
    assert_eq!(
        o(
            "--fastq_mergepairs f --reverse r --fastqout m --log l",
            "log"
        ),
        ("w/log".to_string(), "log-report".to_string())
    );
    // elsewhere --log is a log of the run, not compared
    assert_eq!(
        outputs("--fastx_revcomp in --log l --fastqout o"),
        HashMap::from([(
            "fastqout".to_string(),
            ("w/tmp1.revcomp".to_string(), "bytes".to_string())
        )])
    );
}

#[test]
fn mergepairs_outputs() {
    assert_eq!(
        render(
            "--fastq_mergepairs in --reverse rev --fastq_minovlen 5 --fastqout m \
             --fastqout_notmerged_fwd n --eetabbedout e --log l"
        ),
        "mergpa -i in --reverse rev --min-overlap 5 --not-merged w/tmp1.notmerged_fwd \
         --ee-table w/eetabbedout --report w/log > w/tmp2.merged"
    );
}

#[test]
fn args() {
    let a = Args::parse(&["vsearch", "--fastx_filter=in", "--fastq_maxee=1", "-quiet"]).unwrap();
    assert_eq!(
        (a.command.as_str(), a.input.as_deref()),
        ("fastx_filter", Some("in"))
    );
    assert_eq!(
        a.opts,
        HashMap::from([
            ("fastq_maxee".to_string(), Value::Text("1".into())),
            ("quiet".to_string(), Value::Flag)
        ])
    );
    let a = Args::parse(&["--shuffle", "in", "--randseed", "-5"]).unwrap();
    assert_eq!(a.opts["randseed"], Value::Text("-5".into()));
}

#[test]
fn every_command() {
    let minimal: HashMap<&str, &str> = HashMap::from([
        ("usearch_global", "--db d --id 0.9 --uc u"),
        ("search_exact", "--db d --uc u"),
        ("orient", "--db d --tabbedout t"),
        ("sintax", "--db d --tabbedout t"),
        ("uchime_ref", "--db d --nonchimeras n"),
        ("allpairs_global", "--id 0.9 --uc u"),
        ("cluster_fast", "--id 0.9 --centroids c"),
        ("cluster_size", "--id 0.9 --uc u"),
        ("cluster_smallmem", "--usersort --id 0.9 --uc u"),
        ("cluster_unoise", "--centroids c"),
        ("cut", "--cut_pattern G^AATT_C --fastaout o"),
        ("fastq_join", "--reverse r --fastqout o"),
        ("fastq_mergepairs", "--reverse r --fastqout o"),
        ("fastx_syncpairs", "--reverse r --fastqout o"),
        ("fastx_getsubseq", "--label a --fastaout o"),
        ("fastx_getseq", "--label a --fastaout o"),
        ("fastx_getseqs", "--label_word a --fastaout o"),
        (
            "fastx_subsample",
            "--sample_pct 50 --randseed 1 --fastaout o",
        ),
    ]);
    for cmd in COMMANDS
        .iter()
        .filter(|c| !["help", "version", "h", "v"].contains(c))
    {
        let extra = minimal.get(cmd).copied().unwrap_or("--output o");
        let first = translate(&words(&format!("--{cmd} in {extra}")), W);
        let p = match first {
            Ok(p) => p,
            // commands whose outputs have other names
            Err(_) => [
                "--fastaout o",
                "--fastqout o",
                "--nonchimeras n",
                "--tabbedout t",
                "--log l",
            ]
            .iter()
            .filter_map(|out| translate(&words(&format!("--{cmd} in {out}")), W).ok())
            .find(|p| !p.outputs.is_empty())
            .unwrap_or_else(|| panic!("--{cmd} does not translate")),
        };
        assert!(!p.steps.is_empty(), "{cmd}");
        assert!(!p.outputs.is_empty(), "{cmd}");
    }
}

// compressed inputs go through gzip -dc or bzip2 -dc, known by their magic
// numbers as vsearch knows them

struct Files {
    dir: std::path::PathBuf,
}

impl Files {
    fn new(name: &str) -> Files {
        let dir =
            std::env::temp_dir().join(format!("rscompare-test-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Files { dir }
    }

    /// A file compressed by `tool` (gzip, bzip2), or plain.
    fn write(&self, name: &str, data: &str, tool: Option<&str>) -> String {
        let path = self.dir.join(name);
        match tool {
            None => std::fs::write(&path, data).unwrap(),
            Some(tool) => {
                let plain = self.dir.join(format!("{name}.plain"));
                std::fs::write(&plain, data).unwrap();
                let out = std::process::Command::new(tool)
                    .arg("-c")
                    .arg(&plain)
                    .output()
                    .unwrap();
                std::fs::write(&path, out.stdout).unwrap();
            }
        }
        path.to_string_lossy().into_owned()
    }
}

const FQ: &str = "@a\nACGT\n+\nIIII\n";
const FA: &str = ">d\nACGTACGT\n";

#[test]
fn gzip() {
    let f = Files::new("gzip");
    let q = f.write("q.fq.gz", FQ, Some("gzip"));
    let p = pipe(&format!("--fastx_revcomp {q} --fastqout o"));
    assert_eq!(
        p.render(),
        format!("gzip -dc {q} | revcomp > w/tmp2.revcomp")
    );
    assert!(p.steps[0].system);
    assert_eq!(
        p.steps[0].argv(Path::new("bin")),
        ["gzip", "-dc", q.as_str()]
    );
    assert_eq!(p.steps[1].argv(Path::new("bin")), ["bin/revcomp"]);
}

#[test]
fn format_inside() {
    // FASTA inside the gzip file: no conv at --fasta_width 0
    let f = Files::new("inside");
    let a = f.write("a.fa.gz", FA, Some("gzip"));
    assert_eq!(
        render(&format!("--fastx_revcomp {a} --fasta_width 0 --fastaout o")),
        format!("gzip -dc {a} | revcomp > w/tmp2.revcomp")
    );
}

#[test]
fn database_and_pairs() {
    let f = Files::new("db");
    let q = f.write("q.fq.gz", FQ, Some("gzip"));
    let db = f.write("db.fa.bz2", FA, Some("bzip2"));
    let r = render(&format!(
        "--usearch_global {q} --db {db} --id 1 --qmask soft --dbmask soft --uc u"
    ));
    has(&r, &format!("gzip -dc {q} > w/tmp1.fastq"));
    has(
        &r,
        &format!("bzip2 -dc {db} | sqfilt --min-len 32 --max-len 50000 > w/tmp3.db"),
    );
    has(
        &r,
        "dbsrch --db w/tmp3.db --id 1 --maxhits 1 < w/tmp1.fastq",
    );
    let r = render(&format!(
        "--fastq_mergepairs {q} --reverse {q} --fastqout m"
    ));
    has(&r, "mergpa --reverse w/tmp2.fastq < w/tmp1.fastq");
}

#[test]
fn decompress_flag() {
    // vsearch reads an uncompressed file as it is with --gzip_decompress
    let f = Files::new("flag");
    let plain = f.write("q.fq", FQ, None);
    let q = f.write("q.fq.gz", FQ, Some("gzip"));
    assert_eq!(
        render(&format!(
            "--fastx_revcomp {plain} --gzip_decompress --fastqout o"
        )),
        format!("revcomp -i {plain} > w/tmp1.revcomp")
    );
    assert!(
        render(&format!(
            "--fastx_revcomp {q} --gzip_decompress --fastqout o"
        ))
        .starts_with(&format!("gzip -dc {q} |"))
    );
}

// errors and caveats

#[test]
fn unknown() {
    bad("--fastx_filter in --no_such_option 3", "--no_such_option");
    bad("--fastx_filter in extra", "unexpected argument 'extra'");
    bad("--fastq_maxee 1", "no vsearch command");
    bad("--help", "no rsearch translation");
}

#[test]
fn untranslatable() {
    bad(
        "--fastx_filter in --gzip_decompress --bzip2_decompress --fastqout o",
        "--bzip2_decompress",
    );
    bad("--cluster_fast in --id 0.9 --alnout a", "--alnout with");
    bad(
        "--cluster_fast in --id 0.9 --clusters c",
        "--fasta_width 80",
    );
    bad(
        "--fastx_filter in --fastq_minlen 0 --fastqout o",
        "empty read",
    );
    bad(
        "--fastx_uniques in --minuniquesize 2 --fastaout o",
        "without --sizeout",
    );
    bad(
        "--fastx_filter in --relabel R --fastqout o --fastqout_discarded d",
        "--relabel",
    );
    bad(
        "--fastq_mergepairs in --reverse r --fastq_maxee 1 --fastqout_notmerged_fwd n",
        "not written as not merged",
    );
    bad("--usearch_global q --db d --uc u", "needs --id");
    bad("--usearch_global q --db d --id 1 --userout u", "userfields");
    bad(
        "--usearch_global q --db d --id 1 --fastapairs p",
        "--fasta_width 0",
    );
    bad(
        "--orient q --db d --qmask none --fastaout o",
        "--qmask none",
    );
    bad(
        "--sortbysize in --minseqlength 0 --minsize 2 --output o",
        "empty sequences",
    );
    bad(
        "--fastx_subsample in --sizein --sample_pct 5 --fastaout o",
        "--by-size",
    );
    bad(
        "--uchime_denovo in --fasta_score --nonchimeras n",
        "--fasta_score",
    );
    bad(
        "--fastx_revcomp in --fastq_ascii 64 --fastqout o",
        "--fastq_ascii",
    );
}

#[test]
fn caveat_lists() {
    let c = caveats("--fastx_uniques in --minseqlength 0 --fastaout o");
    assert!(
        c.iter()
            .any(|x| x.contains("--maxseqlength 50000 not applied"))
    );
    assert!(c.iter().any(|x| x.contains("--notrunclabels")));
    assert!(
        caveats("--fastx_uniques in --minseqlength 1 --sizein --notrunclabels --fastaout o")
            .is_empty()
    );
    // vsearch cuts these labels at spaces too, the abundance after one
    // included
    for cmd in [
        "--rereplicate in --output o",
        "--fastx_getseq in --label a --fastaout o",
        "--fastx_getseqs in --label_word a --fastaout o",
        "--fastx_getsubseq in --label a --subseq_start 2 --fastaout o",
    ] {
        assert!(
            caveats(cmd).iter().any(|x| x.contains("--notrunclabels")),
            "{cmd}"
        );
    }
    // vsearch dusts each strand of a query on its own; lcmask masks the plus
    // strand, which dbsrch turns round
    let each = |c: &str| caveats(c).iter().any(|x| x.contains("each strand"));
    assert!(each(
        "--usearch_global q --db d --id 0.9 --strand both --uc u"
    ));
    assert!(!each("--usearch_global q --db d --id 0.9 --uc u"));
    assert!(!each(
        "--usearch_global q --db d --id 0.9 --strand both --qmask soft --uc u"
    ));
    assert!(
        !caveats("--fastx_syncpairs in --reverse r --fastqout o")
            .iter()
            .any(|x| x.contains("--notrunclabels"))
    );
    assert!(
        caveats("--sintax q --db d --tabbedout t")
            .iter()
            .any(|x| x.contains("threads"))
    );
    assert!(
        !caveats("--sintax q --db d --tabbedout t --threads 1 --randseed 1")
            .iter()
            .any(|x| x.contains("threads"))
    );
    assert!(
        caveats("--usearch_global q --db d --id 1 --band 3 --uc u")
            .iter()
            .any(|x| x.contains("ignores"))
    );
}

// threads

/// the arguments of the step running `tool`.
fn step_args(p: &Pipeline, tool: &str) -> Vec<String> {
    p.steps
        .iter()
        .find(|s| s.tool == tool)
        .unwrap_or_else(|| panic!("no {tool} in\n{}", p.render()))
        .args
        .clone()
}

#[test]
fn threads_go_to_the_threaded_tools() {
    let search = "--usearch_global q --db d --id 0.9 --uc out --notrunclabels";
    for (given, expected) in [("1", "1"), ("4", "4"), ("0", "0")] {
        let p = pipe(&format!("{search} --threads {given}"));
        let args = step_args(&p, "dbsrch");
        assert_eq!(&args[args.len() - 2..], ["-t", expected], "{}", p.render());
        // the threaded tools of the pipeline take them, the others none
        for s in &p.steps {
            let threaded = THREADED_TOOLS.contains(&s.tool.as_str());
            assert_eq!(
                s.args.contains(&"-t".to_string()),
                threaded,
                "{}",
                p.render()
            );
        }
    }
    // without --threads, both take one per core
    let p = pipe(search);
    assert!(!step_args(&p, "dbsrch").contains(&"-t".to_string()));
}

#[test]
fn threads_are_set_for_rsearch_alone() {
    let mut p = pipe("--usearch_global q --db d --id 0.9 --uc out --threads 1 --notrunclabels");
    p.set_threads(6);
    let args = step_args(&p, "dbsrch");
    assert_eq!(args.iter().filter(|a| *a == "-t").count(), 1);
    assert_eq!(&args[args.len() - 2..], ["-t", "6"]);
    let mut p = pipe("--usearch_global q --db d --id 0.9 --uc out --notrunclabels");
    p.set_threads(3);
    assert_eq!(
        &step_args(&p, "dbsrch")[step_args(&p, "dbsrch").len() - 2..],
        ["-t", "3"]
    );
}

#[test]
fn threads_reorder_queries_of_some_commands() {
    let order = "vsearch writes each query's records in the order its threads finish";
    let search = "--usearch_global q --db d --id 0.9 --uc out --notrunclabels";
    assert!(
        caveats(&format!("{search} --threads 4"))
            .iter()
            .any(|c| c.contains(order))
    );
    assert!(caveats(search).iter().any(|c| c.contains(order)));
    assert!(
        !caveats(&format!("{search} --threads 1"))
            .iter()
            .any(|c| c.contains(order))
    );
    // clustering gives the same order on any number of threads
    let cluster = "--cluster_fast q --id 0.9 --centroids out --notrunclabels --threads 4";
    assert!(!caveats(cluster).iter().any(|c| c.contains("threads")));
}
