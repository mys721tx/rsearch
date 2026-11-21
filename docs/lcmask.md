# lcmask

Mask low-complexity regions with DUST.

## Synopsis

```
lcmask [-i FILE] [-o FILE] [--hard [--no-dust] | --clear]
```

## Description

lcmask is the only tool that masks. By default it upper-cases each
sequence and lower-cases each low-complexity region DUST finds (soft
masking). With `--hard`, each region becomes N instead and the rest keeps
its case. With `--hard --no-dust`, only the residues already in lower case
become N. `--clear` upper-cases every residue, unmasking it. Headers and
FASTQ scores are written unchanged, in the input format, FASTA unwrapped.

The tools that index or align (`dbsrch`, `clust`, `unoise`, `uchime`,
`mchime`, `sintax`, `dbornt`, `mkudb`) take lower case as masked. lcmask
before them chooses their masking, and `lcmask --clear` turns it off.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `--hard`: Mask with N instead of lower case.
- `--no-dust`: Skip DUST and turn only the residues already in lower case
  into N. Needs `--hard`.
- `--clear`: Upper-case every residue, clearing any masking. Not with
  `--hard` or `--no-dust`.

## vsearch

lcmask is `--fastx_mask` and `--maskfasta`, and the `--qmask`, `--dbmask`
and `--hardmask` options of every command, as runs before the tool that
reads the sequences. `--qmask dust` (the default) is `lcmask`, `--hardmask`
`lcmask --hard`, `--qmask soft --hardmask` `lcmask --hard --no-dust`, and
`--qmask none` `lcmask --clear`. `--qmask soft` is no lcmask at all. The
`--min_unmasked_pct` and `--max_unmasked_pct` filters of `--fastx_mask` are
[sqfilt](sqfilt.md) after it. Like the other streaming tools, lcmask writes
the records before an invalid abundance, where `--fastx_mask` reads the
whole file first and writes none.

## Examples

```console
$ cat repeats.fa
>m1
ACGTTGCAAGGTCACACACACACACACACACAGGTACCAT
>m2
ACGTTGCAAGGTAGGTCAGTTCAGGACCTTAGCATACCAT
$ lcmask -i repeats.fa
>m1
ACGTTGCAAGGTcacacacacacacacacacaGGTACCAT
>m2
ACGTTGCAAGGTAGGTCAGTTCAGGACCTTAGCATACCAT
$ lcmask --hard -i repeats.fa | head -2
>m1
ACGTTGCAAGGTNNNNNNNNNNNNNNNNNNNNGGTACCAT
$ lcmask -i repeats.fa | lcmask --clear | head -2
>m1
ACGTTGCAAGGTCACACACACACACACACACAGGTACCAT
```

vsearch's `--fastx_mask repeats.fa --min_unmasked_pct 80 --fastaout out.fa`:

```console
$ lcmask -i repeats.fa | sqfilt --min-upper-pct 80
>m2
ACGTTGCAAGGTAGGTCAGTTCAGGACCTTAGCATACCAT
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[sqfilt](sqfilt.md), [dbsrch](dbsrch.md), [clust](clust.md),
[mkudb](mkudb.md)
