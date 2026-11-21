# hitfmt

Write hits as alignments, BLAST tabular, user fields, SAM, FASTA pairs or
aligned segments.

## Synopsis

```
hitfmt HITS -i QUERIES --db DATABASE [-o FILE]
       (--alignments | --blast6 | --userout FIELDS | --sam | --fasta-pairs |
        --qsegs | --tsegs)
       [--sam-header] [--output-no-hits] [--rowlen N] [--iddef 0..4]
       [--n-mismatch | --exact] [--match S] [--mismatch S]
       [--gapopen P] [--gapext P]
```

## Description

hitfmt rebuilds vsearch's hit outputs from UC hits, the queries and the
database they were found in. The queries (`-i`) are matched to the hits in
order, by header; the targets are found in the database (`--db`, FASTA,
FASTQ or UDB) by their number. Each alignment is rebuilt from its CIGAR and
strand. Its scores and counts are recomputed with the scores and identity
definition given. These must be those the search used, since UC records do
not hold them.

It writes one format a run:

- `--alignments`: vsearch's `--alnout`, rows of `--rowlen` columns;
- `--blast6`: BLAST tabular, twelve fields;
- `--userout FIELDS`: vsearch's user fields, joined by `+` (query, target,
  evalue, id, pctpv, pctgaps, pairs, gaps, qlo, qhi, tlo, thi, pv, ql, tl,
  qs, ts, alnlen, opens, exts, raw, bits, aln, caln, qstrand, tstrand,
  qrow, trow, qframe, tframe, mism, ids, qcov, tcov, id0 to id4, qilo,
  qihi, tilo, tihi, diffs, mid, qseq, tseq, qrowdots, trowdots, qlor, qhir,
  tlor, thir);
- `--sam`: SAM, with `--sam-header` a header first, whose `@PG` line names
  rsearch;
- `--fasta-pairs`: the aligned query and target of each hit, FASTA;
- `--qsegs`, `--tsegs`: the aligned segments of the queries or targets,
  FASTA.

FASTA outputs are unwrapped. Queries without hits are written only with
`--output-no-hits`, in the formats that have a way of writing them.

## Options

- `HITS`: UC file of hits. Required.
- `-i`, `--in FILE`: FASTA or FASTQ file of the queries. Required.
- `--db FILE`: Database: FASTA, FASTQ or UDB file. Required.
- `-o`, `--out FILE`: Output file. Default: standard output.
- `--alignments`: Write alignments, as vsearch's `--alnout`.
- `--blast6`: Write BLAST tabular, as vsearch's `--blast6out`.
- `--userout FIELDS`: Write these user fields, separated by `+`, as
  vsearch's `--userout` with `--userfields`.
- `--sam`: Write SAM, as vsearch's `--samout`.
- `--fasta-pairs`: Write the aligned rows of each hit, as vsearch's
  `--fastapairs`.
- `--qsegs`: Write the aligned segments of the queries, as vsearch's
  `--qsegout`.
- `--tsegs`: Write the aligned segments of the targets, as vsearch's
  `--tsegout`.
- `--sam-header`: Begin SAM with a header.
- `--output-no-hits`: Also write the queries without hits, as vsearch's
  `--output_no_hits`.
- `--rowlen N`: Columns of an alignment row, 0 for one row. Default: 64.
- `--iddef 0..4`: Identity definition, as the search's. Default: 2.
- `--n-mismatch`: Count N as a mismatch, as the search did.
- `--exact`: The hits are those of [exsrch](exsrch.md). Score every
  position a match, ambiguous symbols included. Not with `--n-mismatch`.
- `--match SCORE`: Score of a match, as the search's. Default: 2.
- `--mismatch SCORE`: Score of a mismatch, as the search's. Default: -4.
- `--gapopen PENALTIES`: Gap opening penalties, as the search's.
  Default: 20I/2E.
- `--gapext PENALTIES`: Gap extension penalties, as the search's.
  Default: 2I/1E.

## vsearch

hitfmt stands for the `--alnout`, `--blast6out`, `--userout`, `--samout`,
`--fastapairs`, `--qsegout` and `--tsegout` outputs of vsearch's search
commands (`--usearch_global`, `--search_exact`, `--allpairs_global`). It
writes no opening lines before the alignments of `--alignments`. vsearch's
`--search_exact` scores an exact hit as its length in matches, an N facing
an N included, where its aligner scores those 0. A UC hit does not say
which search made it, and `--exact` tells hitfmt. Where an opening penalty
is below its extension, vsearch's aligner scores gaps at the right ends
otherwise than the CIGAR gives. The `raw` field can differ there. Two
adjacent queries with one header, both with hits, read as one.

## Examples

```console
$ dbsrch -i queries.fa --db db.fa --id 0.97 --strand both > hits.uc
$ hitfmt hits.uc -i queries.fa --db db.fa --userout query+id+mism+qstrand
q1;sample=S1;size=3	99.4	1	+
q2;sample=S1;size=1	100.0	0	+
q3;sample=S2;size=2	100.0	0	+
q4;sample=S2;size=1	100.0	0	-
$ hitfmt hits.uc -i queries.fa --db db.fa --alignments --rowlen 40 | head -16

Query >q1;sample=S1;size=3
 %Id   TLen  Target
 99%    160  ref1;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis

 Query 160nt >q1;sample=S1;size=3
Target 160nt >ref1;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis

Qry   1 + GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCT 40
          ||||||||||||||||||||||||||||||||||||||||
Tgt   1 + GCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCT 40

Qry  41 + TATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCA 80
           |||||||||||||||||||||||||||||||||||||||
Tgt  41 + GATTTACCCACTCTGCCAAACTCCAGCGCGGTCAGTTCCA 80

$ exsrch -i queries.fa --db db.fa > exact.uc
$ hitfmt exact.uc -i queries.fa --db db.fa --exact --blast6 | cut -f 1,3-12
q2;sample=S1;size=1	100.0	160	0	0	1	160	1	160	-1	0
q3;sample=S2;size=2	100.0	160	0	0	1	160	1	160	-1	0
```

## Exit status

0 on success, 1 for invalid input (including hits that do not match the
queries or the database given), 2 for invalid arguments.

## See also

[dbsrch](dbsrch.md), [exsrch](exsrch.md), [allpr](allpr.md),
[otutab](otutab.md), [lcatax](lcatax.md)
