# eestat

Report quality scores and expected errors at each position of FASTQ reads.

## Synopsis

```
eestat [-i FILE] [-o FILE] [--ascii 33|64] [--qmin Q] [--qmax Q]
```

## Description

eestat reads FASTQ and writes vsearch's `--fastq_eestats` table, a line for
each position. A line gives the number and share of the reads reaching the
position. It then gives the minimum, lower quartile, median, mean, upper
quartile and maximum of the quality scores at it, of their error
probabilities, and of the expected errors of the reads up to it. Negative
scores count as zero.

## Options

- `-i`, `--in FILE`: Input FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the table. Default: standard output.
- `--ascii 33|64`: ASCII offset of the quality scores. Default: 33.
- `--qmin Q`: Lowest score accepted in the input. Default: 0.
- `--qmax Q`: Highest score accepted in the input. Default: 126 minus
  `--ascii`.

## vsearch

eestat is `--fastq_eestats`, which writes the table to `--output`.

## Examples

```console
$ eestat -i reads.fq | head -3 | cut -f 1-9
Pos	Recs	PctRecs	Min_Q	Low_Q	Med_Q	Mean_Q	Hi_Q	Max_Q
1	4	100.0	30.0	30.0	38.0	37.0	40.0	40.0
2	4	100.0	28.0	28.0	38.0	36.2	39.0	40.0
```

## Exit status

0 on success, 1 for invalid input (including FASTA, or a score out of
range), 2 for invalid arguments.

## See also

[eegrid](eegrid.md), [fqstat](fqstat.md), [sqtrim](sqtrim.md)
