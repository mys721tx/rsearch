# sqtrim

Shorten reads: strip their ends and truncate them.

## Synopsis

```
sqtrim [-i FILE] [-o FILE] [--strip-left N] [--strip-right N] [--len N]
       [--qual Q] [--max-ee EE] [--max-ee-rate RATE]
       [--ascii 33|64] [--qmin Q] [--qmax Q]
```

## Description

sqtrim shortens each read in four steps, in this order:

1. `--strip-left` and `--strip-right` remove bases from the 5' and 3' ends.
2. `--len` truncates the read to a length; shorter reads are kept whole.
3. `--qual` truncates before the first base with a score at or below it.
4. `--max-ee` and `--max-ee-rate` truncate before the expected errors of
   the read, or the expected errors per base, would exceed them.

Reads are never discarded, even when emptied. Discarding is
[sqfilt](sqfilt.md)'s job. Its default `--min-len 1` discards the empty
reads that vsearch writes to its discarded output. Expected errors are the
sum of 10^(-Q/10) over the bases. Records are written in the input format,
FASTA unwrapped. The options that read scores are refused for FASTA input.

## Options

- `-i`, `--in FILE`: Input FASTQ or FASTA file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `--strip-left N`: Remove N bases from the 5' end.
- `--strip-right N`: Remove N bases from the 3' end.
- `--len N`: Truncate to N bases, at least 1; shorter reads are kept.
- `--qual Q`: Truncate before the first base with a score of Q or less,
  0 to 93.
- `--max-ee EE`: Truncate before the expected errors exceed EE.
- `--max-ee-rate RATE`: Truncate before the expected errors per base exceed
  RATE.
- `--ascii 33|64`: ASCII offset of the quality scores. Default: 33.
- `--qmin Q`: Lowest score accepted in the input. Default: 0.
- `--qmax Q`: Highest score accepted in the input. Default: 126 minus
  `--ascii`.

## vsearch

sqtrim is the trimming half of `--fastq_filter` and `--fastx_filter`:
`--fastq_stripleft`, `--fastq_stripright`, `--fastq_trunclen_keep`
(`--len`), `--fastq_truncqual`, `--fastq_truncee` and `--fastq_truncee_rate`.
vsearch's `--fastq_trunclen N` truncates and then discards shorter reads,
which is `sqtrim --len N | sqfilt --min-len N`. vsearch filters after it
trims, and so does the pipeline, since [sqfilt](sqfilt.md) sees sqtrim's
output.

## Examples

```console
$ cat reads.fq
@r1
ACGTACGTAGCTAGCTAGGATCGATCGA
+
IHGFEDCBA@?>=<:9876543210/.-
@r2
TTGACCGTAGCAAGCTTAGCNAGCTAAT
+
GGGFFFEEEDDDCCBBBAAA@@@???>>
@r3
GGCATCGATCGA
+
?=;97531/-+)
@r4
ACGTACGTAGCTAGCTAGGATCGATCGA
+
IIIIIIHHHHHHGGGGGFFFFFFEEEEE
$ sqtrim --strip-left 2 --len 10 -i reads.fq
@r1
GTACGTAGCT
+
GFEDCBA@?>
@r2
GACCGTAGCA
+
GFFFEEEDDD
@r3
CATCGATCGA
+
;97531/-+)
@r4
GTACGTAGCT
+
IIIIHHHHHH
$ sqtrim --qual 15 -i reads.fq | grep -A1 '^@r3'
@r3
GGCATCGA
$ sqtrim --max-ee 0.05 -i reads.fq | awk "NR % 4 == 2"
ACGTACGTAGCTAGCTAGGA
TTGACCGTAGCAAGCTTAGCNAGCTAAT
GGCATCG
ACGTACGTAGCTAGCTAGGATCGATCGA
```

vsearch's `--fastq_filter reads.fq --fastq_trunclen 20 --fastq_maxee 0.5
--fastaout out.fa`:

```console
$ sqtrim --len 20 -i reads.fq | sqfilt --min-len 20 --max-ee 0.5 | conv --fasta
>r1
ACGTACGTAGCTAGCTAGGA
>r2
TTGACCGTAGCAAGCTTAGC
>r4
ACGTACGTAGCTAGCTAGGA
```

## Exit status

0 on success, 1 for invalid input (such as a score outside `--qmin` to
`--qmax`), 2 for invalid arguments:

```console
$ sqtrim --qual 20 -i seqs.fa; echo "exit $?"
error: --qual needs quality scores, which FASTA input lacks

Usage: sqtrim [OPTIONS]

For more information, try '--help'.
exit 2
```

## See also

[sqfilt](sqfilt.md), [conv](conv.md), [mergpa](mergpa.md),
[subseq](subseq.md)
