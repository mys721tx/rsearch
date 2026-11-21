# otutab

Count hits per sample and OTU: classic, mothur or BIOM table.

## Synopsis

```
otutab HITS --db FILE [-o FILE] [--mothur | --biom]
```

## Description

otutab reads UC hits and the database they were found in (`--db`: FASTA,
FASTQ or UDB) and writes a table of counts, an OTU a row and a sample a
column in the classic table:

- each query of the hits counts its abundance (`;size=`, or 1) for the OTU
  of its first hit's target;
- a query's sample is its `sample=` or `barcodelabel=` annotation, or else
  the letters, digits and underscores that begin its label;
- an OTU is the target's `otu=` annotation, or else its label up to the
  first semicolon;
- every database sequence is an OTU, with or without hits, and its `tax=`
  annotation, if any, is its taxonomy (the last label read of an OTU gives
  it).

`--mothur` writes a mothur shared table instead (a sample a row, labelled
`vsearch` as vsearch does). `--biom` writes a BIOM 1.0 table in JSON. Its
`id` is the output file's name (or `-`), its `generated_by` rsearch and its
`date` the time of the run in UTC.

## Options

- `HITS`: UC file of hits. Required.
- `--db FILE`: Database the hits were found in: FASTA, FASTQ or UDB file.
  Required.
- `-o`, `--out FILE`: Output table. Default: standard output.
- `--mothur`: Write a mothur shared table.
- `--biom`: Write a BIOM 1.0 table.

## vsearch

otutab is the `--otutabout`, `--mothur_shared_out` and `--biomout` outputs
of vsearch's search commands, on the hits of [dbsrch](dbsrch.md) or
[exsrch](exsrch.md). vsearch counts a sequence as matched only by accepted
hits, which UC records do not tell from weak ones. When weak hits are
written (`--weak-id`), the taxonomy of an OTU whose sequences carry
different ones can therefore differ.

## Examples

```console
$ dbsrch -i queries.fa --db db.fa --id 0.97 --strand both > hits.uc
$ otutab hits.uc --db db.fa | cut -f 1-3
#OTU ID	S1	S2
ref1	3	0
ref2	0	2
ref3	1	1
$ otutab hits.uc --db db.fa --mothur
label	Group	numOtus	ref1	ref2	ref3
vsearch	S1	3	3	0	1
vsearch	S2	3	0	2	1
```

## Exit status

0 on success, 1 for invalid input (including hits naming a target the
database lacks), 2 for invalid arguments.

## See also

[dbsrch](dbsrch.md), [exsrch](exsrch.md), [lcatax](lcatax.md),
[relab](relab.md)
