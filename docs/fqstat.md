# fqstat

Report the lengths and quality scores of FASTQ reads.

## Synopsis

```
fqstat [-i FILE] [-o FILE] [--ascii 33|64] [--qmin Q] [--qmax Q]
```

## Description

fqstat reads FASTQ and writes vsearch's `--fastq_stats` report, in five
tables and a summary:

1. the distribution of read lengths;
2. the distribution of quality scores, with their error probabilities;
3. for each length L, the share of reads at least that long, and the
   average score, error probability and expected errors of their first L
   bases;
4. for each length L, the number and share of the reads whose first L
   bases have at most 1, 0.5, 0.25 and 0.1 expected errors;
5. for each length L, the share of the reads that are still at least L
   long when truncated before their first base with a score of 5, 10, 15
   or 20 or less;
6. the number of reads, their average length and the number of bases.

## Options

- `-i`, `--in FILE`: Input FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the report. Default: standard output.
- `--ascii 33|64`: ASCII offset of the quality scores. Default: 33.
- `--qmin Q`: Lowest score accepted in the input. Default: 0.
- `--qmax Q`: Highest score accepted in the input. Default: 126 minus
  `--ascii`.

## vsearch

fqstat is `--fastq_stats`, with `--fastq_ascii`, `--fastq_qmin` and
`--fastq_qmax`. vsearch writes the report to `--log`; fqstat writes it to
`-o`.

## Examples

```console
$ fqstat -i reads.fq | head -12

Read length distribution
      L           N      Pct   AccPct
-------  ----------  -------  -------
>=   28           3    75.0%    75.0%
     12           1    25.0%   100.0%

Q score distribution
ASCII    Q       Pe           N      Pct   AccPct
-----  ---  -------  ----------  -------  -------
    I   40  0.00010           7     7.3%     7.3%
    H   39  0.00013           7     7.3%    14.6%
```

## Exit status

0 on success, 1 for invalid input (including FASTA, or a score out of
range), 2 for invalid arguments.

## See also

[fqchar](fqchar.md), [eestat](eestat.md), [eegrid](eegrid.md),
[sqtrim](sqtrim.md), [sqfilt](sqfilt.md)
