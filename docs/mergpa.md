# mergpa

Merge overlapping read pairs.

## Synopsis

```
mergpa [-i FILE] --reverse FILE [-o FILE] [--not-merged FILE]
       [--rev-not-merged FILE] [--report FILE] [--ee-table FILE]
       [--min-overlap N] [--max-diffs N] [--max-diff-pct PCT] [--stagger]
       [--ascii 33|64] [--qmin Q] [--qmax Q] [--qmin-out Q] [--qmax-out Q]
```

## Description

mergpa merges each forward read with the reverse complement of the read at
the same position in `--reverse`, where they overlap. Overlaps are found by
shared 5-mers and scored by ungapped alignment. A pair merges when one
overlap alone scores high enough, is at least `--min-overlap` long and has
at most `--max-diffs` differences and `--max-diff-pct` percent of them. The
merged bases get posterior scores, limited to `--qmin-out` and
`--qmax-out`. The merged read keeps the forward header, whole.

Pairs not merged go to `--not-merged` and `--rev-not-merged`, upper-cased,
each N scored 0. `--report` writes the counts of pairs merged and not, and
why. `--ee-table` writes, for each merged pair, the expected errors of the
overlap in each read and the differences in it. Both inputs are FASTQ with
as many reads each. They are read in chunks of 500 pairs that are checked
before they are merged. Errors thus come in vsearch's order.

A pair with an empty read is never merged, nor one that would merge into
more than a million bases. Nothing is written to standard error.

## Options

- `-i`, `--in FILE`: Input FASTQ file of the forward reads. Default:
  standard input.
- `-o`, `--out FILE`: Output FASTQ file of the merged reads. Default:
  standard output.
- `--reverse FILE`: Input FASTQ file of the reverse reads. Required.
- `--not-merged FILE`: Output file of the forward reads of the pairs not
  merged.
- `--rev-not-merged FILE`: Output file of the reverse reads of the pairs
  not merged.
- `--report FILE`: Output file of the counts of pairs merged and not, and
  why.
- `--ee-table FILE`: Output file of the expected errors of the overlap in
  each read and the differences of each merged pair.
- `--min-overlap N`: Shortest overlap, at least 5. Default: 10.
- `--max-diffs N`: Most differences in the overlap. Default: 10.
- `--max-diff-pct PCT`: Highest percentage of differences in the overlap.
  Default: 100.
- `--stagger`: Also merge pairs whose reverse read runs past the start of
  the forward read; the overhang is cut.
- `--ascii 33|64`: ASCII offset of the quality scores, read and written.
  Default: 33.
- `--qmin Q`: Lowest input score accepted. Default: 0.
- `--qmax Q`: Highest input score accepted. Default: 126 minus `--ascii`.
- `--qmin-out Q`: Lowest score of a merged base. Default: 0.
- `--qmax-out Q`: Highest score of a merged base. Default: 41.

## vsearch

mergpa is `--fastq_mergepairs`, and only merges. `--fastq_minovlen`,
`--fastq_maxdiffs`, `--fastq_maxdiffpct` and `--fastq_allowmergestagger`
are its options, with the score ranges. vsearch's read filters run before
it as [sqtrim](sqtrim.md) and [sqfilt](sqfilt.md) (`--fastq_truncqual`,
`--fastq_minlen`, `--fastq_maxlen`, `--fastq_maxns`). Its merged-read
filters run after it as [sqfilt](sqfilt.md) (`--fastq_maxee`,
`--fastq_minmergelen`, `--fastq_maxmergelen`). vsearch's summary goes to
`--report`, `--eetabbedout` to `--ee-table`.

The reads sqfilt discards before mergpa, and the merged reads it discards
after, are those vsearch writes as not merged. mergpa's not-merged outputs
hold only the pairs that reached it. vsearch scores a base below Q2 as an
error of 0.75 when it applies `--fastq_maxee` to a merged read, where
sqfilt uses 10^(-Q/10). `sqfilt --max-ee` after mergpa can therefore differ
from `--fastq_maxee` on reads with such scores. What vsearch writes before an
error depends on its threads, and is not matched.

## Examples

```console
$ mergpa -i fwd.fq --reverse rev.fq --report report.txt
@p1/1
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCG
+
IIIIHHHGGGFFFFEEEDDDJJJJJJJJJJJJJJJJJJJJAABBBCCCDDDEEEFFFGGG
@p2/1
CTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCACCCCATCGGAC
+
IIIIHHHGGGFFFFEEEDDDJJJJJJJJJAJJJJJJJJJJAABBBCCCDDDEEEFFFGGG
@p3/1
TGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGC
+
IIIIHHHGGGFFFFEEEDDDJJJJJJJJJJJJJJJJJJJJAABBBCCCDDDEEEFFFGGG
$ cat report.txt
         3  Pairs
         3  Merged (100.0%)
         0  Not merged (0.0%)

Statistics of all reads:
     40.00  Mean read length

Statistics of merged reads:
     60.00  Mean fragment length
      0.00  Standard deviation of fragment length
      0.02  Mean expected error in forward sequences
      0.29  Mean expected error in reverse sequences
      0.01  Mean expected error in merged sequences
      0.00  Mean observed errors in merged region of forward sequences
      0.33  Mean observed errors in merged region of reverse sequences
      0.33  Mean observed errors in merged region
```

vsearch's `--fastq_mergepairs fwd.fq --reverse rev.fq --fastq_maxee 0.02
--fastqout merged.fq`:

```console
$ mergpa -i fwd.fq --reverse rev.fq | sqfilt --max-ee 0.02 | grep "^@"
@p1/1
@p2/1
@p3/1
```

## Exit status

0 on success, 1 for invalid input (a malformed record, a score out of
range, files of different lengths), 2 for invalid arguments.

```console
$ mergpa -i fwd.fq --reverse seqs.fa; echo "exit $?"
Error: FASTQ file expected, FASTA file found
exit 1
```

## See also

[joinpa](joinpa.md), [syncpa](syncpa.md), [sqtrim](sqtrim.md),
[sqfilt](sqfilt.md), [relab](relab.md)
