# exsrch

Find exact full-length matches in a database by hashing; write hits.

## Synopsis

```
exsrch --db FILE [-i FILE] [-o FILE] [--strand plus|both] [--maxhits N]
       [--self] [--max-qsize N] [--min-tsize N] [--min-size-ratio RATIO]
       [--max-size-ratio RATIO]
```

## Description

exsrch reads the database (FASTA or FASTQ) whole, then each query in turn.
The hits of a query are the database sequences identical to it, ignoring
case and taking U for T, by target number; with `--strand both`, also
those identical to its reverse complement. Every hit is written in UC
format, as by [dbsrch](dbsrch.md). A query without hits gets an N record.
It does not align and takes no alignment options.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file of queries. Default:
  standard input.
- `-o`, `--out FILE`: Output UC file of hits. Default: standard output.
- `--db FILE`: Database: FASTA or FASTQ file. Required.
- `--strand plus|both`: Search the plus strand, or both. Default: plus.
- `--maxhits N`: Write at most N hits of a query, 1 or more. Default: all.
- `--self`: Skip targets with the query's header.
- `--max-qsize N`: Skip queries of an abundance (`;size=`) above N.
- `--min-tsize N`: Skip targets of an abundance below N.
- `--min-size-ratio RATIO`: Skip targets unless the query's abundance is
  at least RATIO times theirs.
- `--max-size-ratio RATIO`: Skip targets unless the query's abundance is at
  most RATIO times theirs.

## vsearch

exsrch is `--search_exact` with `--uc --uc_allhits`; `--maxhits 1` gives
its plain `--uc`. vsearch reads `--maxhits 0` as no limit, which here is
leaving the option out. Its other acceptance options do nothing or reject
everything for exact matches, or filter by length, which is
[sqfilt](sqfilt.md)'s job. Its database length limit is `sqfilt --max-len
50000` on the database. Its other outputs are [hitfmt](hitfmt.md) (with
`--exact`), [otutab](otutab.md), [lcatax](lcatax.md) and
[getseq](getseq.md) on the hits. With `--strand both`, a query that is its
own reverse complement hits a target on both strands at once. exsrch writes
the plus strand first, where vsearch leaves the order to its C++ library.

## Examples

```console
$ exsrch -i queries.fa --db db.fa --strand both | cut -f 1-5,8-9
N	*	*	*	.	*	q1;sample=S1;size=3
H	2	160	100.0	+	=	q2;sample=S1;size=1
H	1	160	100.0	+	=	q3;sample=S2;size=2
H	2	160	100.0	-	=	q4;sample=S2;size=1
N	*	*	*	.	*	q5;sample=S2;size=1
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[dbsrch](dbsrch.md), [hitfmt](hitfmt.md), [derep](derep.md),
[getseq](getseq.md)
