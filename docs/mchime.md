# mchime

Detect chimeras of two or more parents in long exact sequences.

## Synopsis

```
mchime [-i FILE] [-o FILE] [--chimeras FILE] [--table FILE]
       [--alignments FILE] [--width N] [--parts N] [--max-parents N]
       [--min-segment N] [--diff-pct PCT] [--abskew SKEW]
```

## Description

mchime finds chimeras of up to `--max-parents` parents, as vsearch
`--chimeras_denovo` does, in input sorted by decreasing abundance
(`sortby --size`). It reads every sequence first. Each is checked against
the earlier non-chimeras at least `--abskew` times as abundant. The
non-chimeras become parents in turn. A sequence's `--parts` are searched for
parents (words of 8, the default scores). It is a chimera when the longest
segments matching them, taken in turn, cover it whole with two parents or
more. Each segment must be at least `--min-segment` long with at most
`--diff-pct` percent differences.

The non-chimeras go to `-o` and the chimeras to `--chimeras`, in the input
format and unchanged. `--table` writes vsearch's table of the chimeras, and
`--alignments` their alignments with the parents.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the non-chimeras, in the input format.
  Default: standard output.
- `--chimeras FILE`: Output file of the chimeras.
- `--table FILE`: Output file of vsearch's `--tabbedout` table of the
  chimeras.
- `--alignments FILE`: Output file of vsearch's `--alnout` alignments of
  the chimeras.
- `--width N`: Width of the alignments' rows, 0 for no wrapping.
  Default: 60.
- `--parts N`: Parts a sequence is searched in, 2 to 100. Default: one for
  each hundred bases.
- `--max-parents N`: Most parents of a chimera, 2 to 20. Default: 3.
- `--min-segment N`: Shortest segment of a parent. Default: 10.
- `--diff-pct PCT`: Percentage of differences a segment may hold, 0 to 50.
  Default: 0.
- `--abskew SKEW`: Parents at least SKEW times as abundant. Default: 1.

## vsearch

mchime is `--chimeras_denovo`, with `--chimeras_parts` as `--parts`,
`--chimeras_parents_max` as `--max-parents` (vsearch finds none below 2 and
stops at 20), `--chimeras_length_min` as `--min-segment`,
`--chimeras_diff_pct` as `--diff-pct`, `--tabbedout` as `--table`,
`--alnout` as `--alignments` and `--alignwidth` as `--width`. Its sorting is
[sortby](sortby.md) before it, its relabelling and `--sizeout`
[relab](relab.md) after it.

## Examples

```console
$ sortby --size -i long.fa | mchime --chimeras chimeras.fa --table table.tsv | grep '>'
>l1;size=50
>l2;size=40
>l3;size=30
$ grep '>' chimeras.fa
>l4;size=3
$ cut -f 2-5,18 table.tsv
l4;size=3	l1;size=50	l2;size=40	l3;size=30	Y
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[uchime](uchime.md), [sortby](sortby.md)
