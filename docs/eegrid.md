# eegrid

Report the FASTQ reads kept by length and expected-error cutoffs.

## Synopsis

```
eegrid [-i FILE] [-o FILE] [--min-len N] [--max-len N] [--len-step N]
       [--max-ee EE,...] [--ascii 33|64] [--qmin Q] [--qmax Q]
```

## Description

eegrid reads FASTQ and writes vsearch's `--fastq_eestats2` grid. The length
cutoffs run from `--min-len` by `--len-step` up to the longest read or
`--max-len`. The expected-error cutoffs are those of `--max-ee`. For each
pair of cutoffs, it counts the reads at least that long whose expected
errors up to that length are within the cutoff. These are the reads that
`sqtrim --len L | sqfilt --min-len L --max-ee EE` would keep. A first line
gives the number of reads and their longest and average lengths.

## Options

- `-i`, `--in FILE`: Input FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the grid. Default: standard output.
- `--min-len N`: First length cutoff. Default: 50.
- `--max-len N`: Last length cutoff. Default: the longest read.
- `--len-step N`: Step between length cutoffs. Default: 50.
- `--max-ee EE,...`: Expected-error cutoffs, comma-separated, in the order
  given. Default: 0.5,1,2.
- `--ascii 33|64`: ASCII offset of the quality scores. Default: 33.
- `--qmin Q`: Lowest score accepted in the input. Default: 0.
- `--qmax Q`: Highest score accepted in the input. Default: 126 minus
  `--ascii`.

## vsearch

eegrid is `--fastq_eestats2`. vsearch's `--length_cutoffs 50,*,50` (first,
last or `*`, step) is split into `--min-len`, `--max-len` and `--len-step`,
and `--ee_cutoffs` is `--max-ee`. vsearch writes the grid to `--output` or
`--log`; eegrid writes it to `-o`.

## Examples

```console
$ eegrid --min-len 10 --len-step 10 --max-ee 0.1,0.5 -i reads.fq
4 reads, max len 28, avg 24.0

Length         MaxEE 0.10         MaxEE 0.50
------   ----------------   ----------------
    10          3( 75.0%)          4(100.0%)
    20          3( 75.0%)          3( 75.0%)
```

## Exit status

0 on success, 1 for invalid input (including FASTA, or a score out of
range), 2 for invalid arguments.

## See also

[eestat](eestat.md), [fqstat](fqstat.md), [sqtrim](sqtrim.md),
[sqfilt](sqfilt.md)
