# dbsrch

Search queries against a database; write hits.

## Synopsis

```
dbsrch --db FILE --id SHARE [-i FILE] [-o FILE] [--strand plus|both]
       [--maxaccepts N] [--maxrejects N] [--word-len N]
       [--min-word-matches N] [criteria] [scores] [--maxhits N]
       [--top-hits-only] [-t N]
```

## Description

dbsrch searches each query against a database by global alignment, as
vsearch `--usearch_global` does:

1. The database (FASTA, FASTQ or UDB, recognised by its first bytes) is read
   whole and its words of `--word-len` residues indexed; words holding a
   lower-case (masked) residue are skipped.
2. For each query, in turn, the targets sharing at least
   `--min-word-matches` words with it are ranked by the words they share.
3. They are aligned globally with the query, best first, in batches of
   eight. Each is accepted or rejected by the criteria, until
   `--maxaccepts` are accepted or `--maxrejects` rejected (0 for no
   limit). With
   `--strand both`, the reverse complement of the query is searched too.

The accepted hits and the weak ones (`--weak-id`) are written in UC format,
by decreasing identity, then target. A query without hits gets an N
record. Each H record has ten tab-separated fields: H, the target's number in
the database from 0, the query's length, the identity, the strand, 0, 0,
the CIGAR of the alignment (`=` for a perfect one), the query's header and
the target's. [hitfmt](hitfmt.md), [otutab](otutab.md),
[lcatax](lcatax.md) and [getseq](getseq.md) read these hits.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file of queries. Default:
  standard input.
- `-o`, `--out FILE`: Output UC file of hits. Default: standard output.
- `-t`, `--threads N`: Search on N threads, 0 for one per logical core.
  The hits are the same, in the same order, on any number. Default: 0.
- `--db FILE`: Database: FASTA, FASTQ or UDB file. Required.
- `--strand plus|both`: Search the plus strand of the queries, or both.
  Default: plus.
- `--maxaccepts N`: Stop after N accepted hits on a strand, 0 for no limit.
  Default: 1.
- `--maxrejects N`: Stop after N rejected targets on a strand, 0 for no
  limit. Default: 32.
- `--word-len N`: Length of the words indexed, 3 to 15; a UDB database's
  own length wins, with a warning when they differ. Default: 8.
- `--min-word-matches N`: Words a target must share with the query to be
  aligned. Default: taken from `--word-len` as given, as vsearch's.

Acceptance criteria (`--id` is required):

- `--id SHARE`: Accept hits of at least this identity, 0 to 1.
- `--iddef 0..4`: Identity definition, as vsearch's: 0 matching columns
  over the shorter sequence (CD-HIT), 1 matching columns over all columns,
  2 the same without terminal gaps, 3 one less the mismatches and gap
  openings over the longer sequence (MBL), 4 as 1 (BLAST). Default: 2.
- `--weak-id SHARE`: Also write weak hits of at least this identity, below
  `--id`; they do not count as accepted.
- `--max-diffs N`: Reject hits with more substitutions and interior gap
  columns.
- `--max-subs N`: Reject hits with more substitutions.
- `--max-gaps N`: Reject hits with more interior gaps.
- `--min-cols N`: Reject hits with fewer interior columns.
- `--left-just`: Reject hits with a terminal gap on the left.
- `--right-just`: Reject hits with a terminal gap on the right.
- `--query-cov SHARE`: Reject hits covering less of the query, 0 to 1.
- `--target-cov SHARE`: Reject hits covering less of the target, 0 to 1.
- `--max-id SHARE`: Reject hits of a higher identity, 0 to 1.
- `--mid PCT`: Reject hits whose identity over the letter pairs (gaps left
  out) is lower, in percent.
- `--id-prefix N`: Reject targets whose first N residues differ from the
  query's.
- `--id-suffix N`: Reject targets whose last N residues differ from the
  query's.
- `--max-qt RATIO`: Reject targets unless the query's length is at most
  RATIO times theirs.
- `--min-qt RATIO`: Reject targets unless the query's length is at least
  RATIO times theirs.
- `--min-sl SHARE`: Reject pairs whose shorter length is below this share
  of the longer.
- `--max-sl SHARE`: Reject pairs whose shorter length is above this share
  of the longer.
- `--max-qsize N`: Reject queries of an abundance (`;size=`) above N.
- `--min-tsize N`: Reject targets of an abundance below N.
- `--min-size-ratio RATIO`: Reject targets unless the query's abundance is
  at least RATIO times theirs.
- `--max-size-ratio RATIO`: Reject targets unless the query's abundance is
  at most RATIO times theirs.
- `--self`: Reject targets with the query's header.
- `--self-id`: Reject targets identical to the query.

Scores, as vsearch's:

- `--match SCORE`: Score of a match. Default: 2.
- `--mismatch SCORE`: Score of a mismatch. Default: -4.
- `--gapopen PENALTIES`: Gap opening penalties, as vsearch's: a penalty
  followed by letters for where it applies (I interior, E both ends, L
  left, R right, Q query, T target), several separated by `/`.
  Default: 20I/2E.
- `--gapext PENALTIES`: Gap extension penalties, written as `--gapopen`.
  Default: 2I/1E.
- `--n-mismatch`: Count N as a mismatch rather than as no score.

Hits written:

- `--maxhits N`: Write at most N hits of a query, 1 or more. Default: all.
- `--top-hits-only`: Write only the hits of the top identity.

## vsearch

dbsrch is `--usearch_global` with `--uc --uc_allhits`. It writes every hit,
since the tools that read them need them all. `--maxhits 1` gives vsearch's
plain `--uc`. The other outputs of the command are tools reading the hits:
`--alnout`, `--blast6out`, `--userout`, `--samout`, `--fastapairs`,
`--qsegout` and `--tsegout` are [hitfmt](hitfmt.md); `--otutabout`,
`--mothur_shared_out` and `--biomout` [otutab](otutab.md); `--lcaout`
[lcatax](lcatax.md); `--matched`, `--notmatched`, `--dbmatched` and
`--dbnotmatched` [getseq](getseq.md). The criteria are vsearch's under
names in words (`--maxdiffs` is `--max-diffs`, `--query_cov` `--query-cov`,
`--mintsize` `--min-tsize`, and so on).

vsearch masks both files with DUST by default. It also drops database
sequences shorter than 32 or longer than 50000. In rsearch that is
[lcmask](lcmask.md) on both and [sqfilt](sqfilt.md) on the database before
dbsrch. dbsrch itself only skips the words holding lower case.

With `--strand both`, vsearch reverse-complements each query and then runs
DUST on each strand on its own. DUST does not give the same regions on the
two strands, so a run of one base near the end of a read can be masked on
the minus strand only. lcmask masks the plus strand once, and dbsrch
reverse-complements that, keeping its case. The minus strand then keeps
words that vsearch masks. The ranking of the candidates changes, and with
the default `--maxaccepts 1` another target can be accepted first. dbsrch
keeps the masking it is given, whoever made it, rather than running DUST
again and upper-casing the queries as vsearch does. `--qmask soft` and
`--qmask none` give vsearch's results, and so does `--strand plus`.

vsearch searches on several threads too, and writes each query's hits as its
thread finishes, so the queries come in another order from one run to the
next. dbsrch writes them in input order on any number of threads. That is
the order of vsearch `--threads 1`, with the same lines as vsearch on any
number.

Hits that tie
on vsearch's sort key (the two strands of a target at one identity) are
written by target, plus strand first; vsearch leaves their order to its C++
library.

## Examples

```console
$ dbsrch -i queries.fa --db db.fa --id 0.97 --strand both | cut -f 1-5,8-10
H	0	160	99.4	+	160M	q1;sample=S1;size=3	ref1;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis
H	2	160	100.0	+	=	q2;sample=S1;size=1	ref3;tax=d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli
H	1	160	100.0	+	=	q3;sample=S2;size=2	ref2;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_cereus
H	2	160	100.0	-	=	q4;sample=S2;size=1	ref3;tax=d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli
N	*	*	*	.	*	q5;sample=S2;size=1	*
$ dbsrch -i queries.fa --db db.fa --id 0.97 --weak-id 0.85 --maxaccepts 0 | cut -f 1-4,9-10 | head -3
H	0	160	99.4	q1;sample=S1;size=3	ref1;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis
H	1	160	89.4	q1;sample=S1;size=3	ref2;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_cereus
H	2	160	100.0	q2;sample=S1;size=1	ref3;tax=d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli
```

vsearch's `--usearch_global queries.fa --db db.fa --id 0.97 --blast6out
b.tsv`:

```console
$ lcmask -i db.fa | sqfilt --min-len 32 --max-len 50000 > dbm.fa
$ lcmask -i queries.fa | dbsrch --db dbm.fa --id 0.97 > hits.uc
$ hitfmt hits.uc -i queries.fa --db dbm.fa --blast6 | cut -f 1-4
q1;sample=S1;size=3	ref1;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis	99.4	160
q2;sample=S1;size=1	ref3;tax=d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli	100.0	160
q3;sample=S2;size=2	ref2;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_cereus	100.0	160
```

## Exit status

0 on success, 1 for invalid input (a malformed query or database, a UDB
file that does not read), 2 for invalid arguments, such as a missing
`--id`.

## See also

[exsrch](exsrch.md), [allpr](allpr.md), [hitfmt](hitfmt.md),
[otutab](otutab.md), [lcatax](lcatax.md), [getseq](getseq.md),
[mkudb](mkudb.md), [lcmask](lcmask.md)
