# allpr

Align every pair of sequences in one file; write hits.

## Synopsis

```
allpr (--id SHARE | --accept-all) [-i FILE] [-o FILE] [criteria] [scores]
      [--maxhits N] [--top-hits-only]
```

## Description

allpr aligns each sequence globally with each one after it in the file, as
vsearch `--allpairs_global` does. It uses vsearch's 16-bit aligner, or its
linear-memory one when a score overflows or the sequences are too long for
it. The hits the criteria accept (all of them with `--accept-all`, the
criteria aside) are written in UC format, by decreasing identity, then
target. A sequence without hits gets an N record. The target is named by its
number in the file, from 0. Each hit is written once, under the earlier
sequence as query.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output UC file of hits. Default: standard output.
- `--accept-all`: Accept every hit the criteria do not reject. Not with
  `--id`.
- `--id SHARE`: Accept hits of at least this identity, 0 to 1.
- `--iddef 0..4`: Identity definition, as vsearch's. Default: 2.
- `--weak-id SHARE`: Also write weak hits of at least this identity, below
  `--id`.
- `--max-diffs N`, `--max-subs N`, `--max-gaps N`, `--min-cols N`: Reject
  hits with more differences, substitutions or interior gaps, or fewer
  interior columns.
- `--left-just`, `--right-just`: Reject hits with a terminal gap on that
  side.
- `--query-cov SHARE`, `--target-cov SHARE`: Reject hits covering less of
  the query or target, 0 to 1.
- `--max-id SHARE`: Reject hits of a higher identity, 0 to 1.
- `--mid PCT`: Reject hits whose identity over the letter pairs is lower,
  in percent.
- `--id-prefix N`, `--id-suffix N`: Reject targets whose first or last N
  residues differ from the query's.
- `--max-qt RATIO`, `--min-qt RATIO`: Reject targets unless the query's
  length is at most, or at least, RATIO times theirs.
- `--min-sl SHARE`, `--max-sl SHARE`: Reject pairs whose shorter length is
  below, or above, this share of the longer.
- `--max-qsize N`, `--min-tsize N`: Reject queries of an abundance above N,
  targets of an abundance below N.
- `--min-size-ratio RATIO`, `--max-size-ratio RATIO`: Reject targets unless
  the query's abundance is at least, or at most, RATIO times theirs.
- `--self`: Reject targets with the query's header.
- `--self-id`: Reject targets identical to the query.
- `--match SCORE`: Score of a match. Default: 2.
- `--mismatch SCORE`: Score of a mismatch. Default: -4.
- `--gapopen PENALTIES`: Gap opening penalties, as vsearch's.
  Default: 20I/2E.
- `--gapext PENALTIES`: Gap extension penalties, as vsearch's.
  Default: 2I/1E.
- `--n-mismatch`: Count N as a mismatch.
- `--maxhits N`: Write at most N hits of a query, 1 or more. Default: all.
- `--top-hits-only`: Write only the hits of the top identity.

[dbsrch](dbsrch.md) describes the criteria and scores at more length.

## vsearch

allpr is `--allpairs_global` with `--uc`, which writes every hit. vsearch
lets `--acceptall` win over `--id`; allpr takes one or the other. vsearch's
`--maxaccepts`, `--maxrejects` and word options do nothing there, and are
left out. Like vsearch it keeps the counts of the 16-bit aligner, which
takes every column of an empty query's alignment for a gap (`--iddef 3`
shows it). vsearch's default `--minseqlength 1` is `sqfilt` before it.
Its other outputs are [hitfmt](hitfmt.md) on the hits.

## Examples

```console
$ allpr -i amplicons.fa --id 0.9 | cut -f 1-5,9-10
H	3	160	99.4	+	a1;size=60	a4;size=2
H	4	160	95.0	+	a1;size=60	a5;size=4
H	1	160	90.0	+	a1;size=60	a2;size=40
H	4	160	95.0	+	a2;size=40	a5;size=4
N	*	*	*	.	a3;size=25	*
H	4	160	94.4	+	a4;size=2	a5;size=4
N	*	*	*	.	a5;size=4	*
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[dbsrch](dbsrch.md), [hitfmt](hitfmt.md), [sqfilt](sqfilt.md)
