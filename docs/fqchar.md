# fqchar

Report the residues and score symbols of FASTQ reads, and guess their
encoding.

## Synopsis

```
fqchar [-i FILE] [-o FILE] [--tail N]
```

## Description

fqchar reads FASTQ and writes vsearch's `--fastq_chars` report:

- the number of reads, the lowest and highest score symbols and their
  range, and the encoding and score range they suggest, as vsearch options
  and in words;
- for each residue, upper-cased: its count, its share and its longest run
  (less one, as vsearch counts it), and for N the lowest and highest score
  symbols seen with it;
- for each score symbol: its ASCII code, its share and the number of reads
  ending in a run of at least `--tail` of it.

## Options

- `-i`, `--in FILE`: Input FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the report. Default: standard output.
- `--tail N`: Shortest run of a score symbol at the end of a read that
  counts as a tail. Default: 4.

## vsearch

fqchar is `--fastq_chars`, with `--fastq_tail` as `--tail`. vsearch writes
the report to standard error and `--log`; fqchar writes it to `-o`, without
vsearch's banner and progress lines.

## Examples

```console
$ fqchar -i reads.fq
Read 4 sequences.
Qmin 41, Qmax 73, Range 33
Guess: -fastq_qmin 8 -fastq_qmax 40 -fastq_ascii 33
Guess: Original Sanger format (phred+33)

Letter          N   Freq MaxRun
------ ---------- ------ ------
     A         27  28.1%      1
     C         21  21.9%      1
     G         26  27.1%      1
     N          1   1.0%      0  Q=@
     T         21  21.9%      1

Char  ASCII    Freq       Tails
----  -----  ------  ----------
 ')'     41    1.0%           0
 '+'     43    1.0%           0
 '-'     45    2.1%           0
 '.'     46    1.0%           0
 '/'     47    2.1%           0
 '0'     48    1.0%           0
 '1'     49    2.1%           0
 '2'     50    1.0%           0
 '3'     51    2.1%           0
 '4'     52    1.0%           0
 '5'     53    2.1%           0
 '6'     54    1.0%           0
 '7'     55    2.1%           0
 '8'     56    1.0%           0
 '9'     57    2.1%           0
 ':'     58    1.0%           0
 ';'     59    1.0%           0
 '<'     60    1.0%           0
 '='     61    2.1%           0
 '>'     62    3.1%           0
 '?'     63    5.2%           0
 '@'     64    4.2%           0
 'A'     65    4.2%           0
 'B'     66    4.2%           0
 'C'     67    3.1%           0
 'D'     68    4.2%           0
 'E'     69    9.4%           1
 'F'     70   10.4%           0
 'G'     71    9.4%           0
 'H'     72    7.3%           0
 'I'     73    7.3%           0
```

## Exit status

0 on success, 1 for invalid input (including FASTA), 2 for invalid
arguments.

## See also

[fqstat](fqstat.md), [eestat](eestat.md), [eegrid](eegrid.md),
[conv](conv.md)
