//! Tests of `benches/common/data.rs`, the module through which the
//! benchmarks of both crates read a real dataset instead of their
//! synthetic input.

use std::path::{Path, PathBuf};

#[path = "../../benches/common/data.rs"]
mod data;

use data::*;

#[test]
fn table_paths_are_relative_to_the_dataset() {
    let dir = Path::new("/data/geo");
    let table = parse_table(
        "# role\tfile\nreads\tall_R1.fastq\n\ndb\t/refs/rdp.fa\nderep.reads\tsub/x.fq\n",
        dir,
    );
    assert_eq!(
        table,
        [
            ("reads".to_string(), PathBuf::from("/data/geo/all_R1.fastq")),
            ("db".to_string(), PathBuf::from("/refs/rdp.fa")),
            (
                "derep.reads".to_string(),
                PathBuf::from("/data/geo/sub/x.fq")
            ),
        ]
    );
}

#[test]
fn a_bench_key_overrides_the_role() {
    let d = Dataset {
        id: "x".into(),
        dir: PathBuf::from("/d"),
        files: parse_table("reads\ta.fq\nderep.reads\tb.fq\n", Path::new("/d")),
    };
    assert_eq!(d.lookup("derep", "reads"), Some(Path::new("/d/b.fq")));
    assert_eq!(d.lookup("sqtrim", "reads"), Some(Path::new("/d/a.fq")));
    assert_eq!(d.lookup("sqtrim", "mate"), None);
}

#[test]
fn group_names_carry_the_dataset() {
    assert_eq!(group_for("derep/binary", None), "derep/binary");
    assert_eq!(
        group_for("derep/binary", Some("geo-miseq-16s-v3v4")),
        "derep/binary@geo-miseq-16s-v3v4"
    );
}

#[test]
fn records_of_fasta_and_fastq() {
    assert_eq!(records(""), 0);
    assert_eq!(records(">a\nAC\nGT\n>b\n\n>c\nA\n"), 3);
    assert_eq!(records("@a\nACGT\n+\nIIII\n@b\nA\n+a\nI\n"), 2);
    // a quality line starting with @ is not a record
    assert_eq!(records("@a\nAC\n+\n@I\n"), 1);
}

#[test]
fn parse_joins_fasta_lines_and_keeps_qualities() {
    let fa = parse(">a x;size=2\nAC\nGT\r\n>b\n\n");
    assert_eq!(fa.len(), 2);
    assert_eq!(fa[0].label, "a x;size=2");
    assert_eq!(fa[0].seq, "ACGT");
    assert_eq!(fa[0].qual, None);
    assert_eq!(fa[1].seq, "");
    let fq = parse("@r1\nACGT\n+r1\nIIII\n@r2\nAC\n+\n@@\n");
    assert_eq!(fq.len(), 2);
    assert_eq!(fq[1].label, "r2");
    assert_eq!(fq[1].qual.as_deref(), Some("@@"));
    assert_eq!(to_fastq(&fq), "@r1\nACGT\n+\nIIII\n@r2\nAC\n+\n@@\n");
    assert_eq!(to_fasta(&fq), ">r1\nACGT\n>r2\nAC\n");
}

#[test]
fn phred64_shifts_qualities_only() {
    assert_eq!(phred64("@r\nACGT\n+\n#I5+\n"), "@r\nACGT\n+\nBhTJ\n");
}

#[test]
fn distinct_sequences_in_first_appearance() {
    let rs = parse(">a\nAC\n>b\nGT\n>c\nAC\n>d\nac\n");
    assert_eq!(distinct(&rs), ["AC", "GT", "ac"]);
}

#[test]
fn sizes_sum_the_annotations() {
    assert_eq!(total_size(">a;size=3\nA\n>b;x=1;size=10;\nC\n>c\nG\n"), 13);
}

#[test]
fn expected_errors_and_length_are_appended() {
    let rs = parse("@r;size=2\nAC\n+\n++\n");
    let out = with_ee_length(&rs);
    // Q10 twice: 0.1 + 0.1
    assert_eq!(out[0].label, "r;size=2;ee=0.2000;length=2");
    assert_eq!(out[0].seq, "AC");
}

#[test]
fn shuffle_is_seeded() {
    let mut a: Vec<u32> = (0..50).collect();
    let mut b = a.clone();
    shuffle(&mut a, 7);
    shuffle(&mut b, 7);
    assert_eq!(a, b);
    assert_ne!(a, (0..50).collect::<Vec<_>>());
    a.sort();
    assert_eq!(a, (0..50).collect::<Vec<_>>());
}

#[test]
fn missing_roles_are_reported() {
    let d = Dataset {
        id: "x".into(),
        dir: PathBuf::from("/d"),
        files: parse_table("reads\ta.fq\n", Path::new("/d")),
    };
    assert_eq!(missing_in(None, "mergpa", &["reads", "mate"]), None);
    assert_eq!(missing_in(Some(&d), "fqchar", &["reads"]), None);
    assert_eq!(
        missing_in(Some(&d), "mergpa", &["reads", "mate"]),
        Some("mate".to_string())
    );
}

#[test]
fn exact_databases_hold_distinct_sequences_and_duplicates() {
    let rs = parse(">r1\nAC\n>r2\nGT\n>r3\nAC\n>r4\nTT\n>r5\nCC\n>r6\nGG\n");
    assert_eq!(
        exact_db(&rs),
        ">d0\nAC\n>d1\nGT\n>d2\nTT\n>d3\nCC\n>d4\nGG\n>dup0\nAC\n>dup4\nGG\n"
    );
}

#[test]
fn optional_roles() {
    let dir = std::env::temp_dir().join(format!("bench-data-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("a.sff"), "").unwrap();
    let d = Dataset {
        id: "x".into(),
        dir: dir.clone(),
        files: parse_table("sff\ta.sff\nphred64\tgone.fq\n", &dir),
    };
    assert_eq!(optional_in(None, "conv", "sff"), None);
    assert_eq!(
        optional_in(Some(&d), "conv", "sff"),
        Some(dir.join("a.sff"))
    );
    assert_eq!(optional_in(Some(&d), "conv", "solexa"), None);
    // a listed file that is not there is not used
    assert_eq!(optional_in(Some(&d), "conv", "phred64"), None);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn output_names_carry_the_dataset_only_with_one() {
    let p = PathBuf::from("/t/clmsa-bench.uc");
    assert_eq!(output_for(p.clone(), None), p);
    assert_eq!(
        output_for(p, Some("geo")),
        PathBuf::from("/t/geo-clmsa-bench.uc")
    );
}

#[test]
fn temporary_names_carry_the_dataset() {
    let p = tmp_for("derep", "x.fq", "geo");
    assert!(p.ends_with("derep-geo-x.fq"), "{p:?}");
    assert_eq!(p.parent(), Some(Path::new(env!("CARGO_TARGET_TMPDIR"))));
}
