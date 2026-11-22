# uchime

Detect chimeras with UCHIME, de novo or against a database.

## Synopsis

```
uchime [-i FILE] [-o FILE] [--db FILE] [--chimeras FILE]
       [--borderline FILE] [--table FILE] [--alignments FILE] [--width N]
       [--uchime2] [--abskew SKEW] [--minh SCORE] [--mindiv PCT]
       [--mindiffs N] [--xn WEIGHT] [--dn COUNT] [-t N]
```

## Description

uchime finds chimeras: sequences that are a left part of one parent joined
to the right part of another.

- De novo (without `--db`), as vsearch `--uchime_denovo`, the input is
  taken as sorted by decreasing abundance (`sortby --size`). Every sequence
  is read first. Each is checked against the earlier non-chimeras at least
  `--abskew` times as abundant. The non-chimeras become parents in turn.
- With `--db` (FASTA, FASTQ or UDB), as `--uchime_ref`, the database is
  read whole and each query, in turn, is checked against it.

The four parts of a sequence are searched for parents (words of 8, lower
case masked, the default scores). The two best explain the sequence as a
chimera or not. By default the test is UCHIME's score, which must reach
`--minh` with at least `--mindiv` percent divergence from the closest
parent and `--mindiffs` differences on each side. With `--uchime2`, the
sequence is a chimera when a model of the two parents matches it
perfectly.

The non-chimeras are written to `-o`, the chimeras to `--chimeras`, and the
borderline sequences (scoring at least `--minh` but failing the other
tests) to `--borderline`, in the input format and unchanged. `--table`
writes vsearch's `--uchimeout` table, a line per sequence (score, query,
parents, identities, votes, divergence, Y, N or ?). `--alignments` writes
the alignments of each chimera with its parents.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the non-chimeras, in the input
  format. Default: standard output.
- `--db FILE`: Database of parents, FASTA, FASTQ or UDB, rather than de
  novo.
- `--chimeras FILE`: Output file of the chimeras.
- `--borderline FILE`: Output file of the borderline sequences, scoring at
  least `--minh` but not chimeras.
- `--table FILE`: Output file of vsearch's `--uchimeout` table.
- `--alignments FILE`: Output file of vsearch's `--uchimealns` alignments
  of the chimeras.
- `--width N`: Width of the alignments' rows, 0 for no wrapping.
  Default: 80.
- `--uchime2`: Call chimeras by UCHIME2: a model matching the sequence
  perfectly. Not with `--db`, `--minh`, `--mindiv` or `--mindiffs`.
- `--abskew SKEW`: De novo, parents at least SKEW times as abundant.
  Default: 2. Not with `--db`.
- `--minh SCORE`: Lowest score of a chimera. Default: 0.28.
- `--mindiv PCT`: Lowest divergence of a chimera from its closest parent,
  in percent. Default: 0.8.
- `--mindiffs N`: Fewest differences on each side of a chimera's
  crossover. Default: 3.
- `--xn WEIGHT`: Weight of the no votes. Default: 8.
- `--dn COUNT`: Pseudo-count of the no votes. Default: 1.4.
- `-t`, `--threads N`: Against a database, judge the queries on N threads,
  0 for one per logical core. The outputs are the same, in the same order,
  on any number. De novo, each sequence's parents are the non-chimeras
  before it, and uchime runs on one thread. Default: 0.

## vsearch

uchime is the four UCHIME commands: `--uchime_denovo`, `--uchime2_denovo`
(`--uchime2`), `--uchime3_denovo` (`--uchime2 --abskew 16`) and
`--uchime_ref` (`--db`), with `--uchimeout` as `--table`, `--uchimealns` as
`--alignments` and `--alignwidth` as `--width`. `--uchime2` refuses the
options vsearch ignores with it. The sorting, the length limit
(`--maxseqlength 50000`) and the masking that vsearch applies first are
[sortby](sortby.md), [sqfilt](sqfilt.md) and [lcmask](lcmask.md) before
it. `--sizeout`, `--fasta_score` and relabelling are [relab](relab.md)
after it. `--uchimeout5` is left out.
Against a database, it also reads FASTQ queries, which vsearch does not.

vsearch runs `--uchime_ref` on several threads and writes each query as its
thread finishes, in another order from one run to the next. uchime writes
the queries in input order on any number of threads. That is the order of
vsearch `--threads 1`, with the same records as vsearch on any number.
vsearch runs the de novo commands on one thread, as uchime does.

## Examples

```console
$ sortby --size -i amplicons.fa | uchime --chimeras chimeras.fa --table table.tsv | grep '>'
>a1;size=60
>a2;size=40
>a3;size=25
>a4;size=2
$ grep '>' chimeras.fa
>a5;size=4
$ cut -f 1-5,18 table.tsv
0.0000	a1;size=60	*	*	*	N
0.0000	a2;size=40	*	*	*	N
0.0000	a3;size=25	*	*	*	N
0.5102	a5;size=4	a1;size=60	a2;size=40	a1;size=60	Y
0.0000	a4;size=2	*	*	*	N
```

vsearch's `--uchime3_denovo amplicons.fa --nonchimeras n.fa`, where no parent is
sixteen times as abundant as a5, and `--uchime_ref amplicons.fa --db db.fa
--chimeras c.fa`:

```console
$ sqfilt --max-len 50000 -i amplicons.fa | lcmask | sortby --size | uchime --uchime2 --abskew 16 | grep '>'
>a1;size=60
>a2;size=40
>a3;size=25
>a5;size=4
>a4;size=2
$ uchime --db db.fa -i amplicons.fa --chimeras ref_chimeras.fa > /dev/null
$ grep ">" ref_chimeras.fa
>a5;size=4
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[mchime](mchime.md), [sortby](sortby.md), [unoise](unoise.md),
[sqfilt](sqfilt.md)
