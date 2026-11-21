# redig

Cut sequences at a restriction pattern.

## Synopsis

```
redig [-i FILE] [-o FILE] --pattern PATTERN [--rev-out FILE]
      [--uncut FILE] [--rev-uncut FILE]
```

## Description

redig cuts each sequence where it matches the pattern, a restriction site
in IUPAC codes with `^` where the forward strand is cut and `_` where the
reverse strand is, as `G^AATT_C` for EcoRI. Codes match when they share a
base, and matches may overlap. The fragments of the forward strand go to
`-o` in order, those of the reverse strand to `--rev-out`, from its 3' end,
each with the header of its sequence. Sequences never cut go to `--uncut`,
and their reverse complements to `--rev-uncut`. FASTA in and out, unwrapped.

## Options

- `-i`, `--in FILE`: Input FASTA file. Default: standard input.
- `-o`, `--out FILE`: Output file of the fragments of the forward strand.
  Default: standard output.
- `--pattern PATTERN`: Restriction pattern, with `^` and `_` at the cut
  sites. Required.
- `--rev-out FILE`: Output file of the fragments of the reverse strand.
- `--uncut FILE`: Output file of the sequences never cut.
- `--rev-uncut FILE`: Output file of the reverse complements of the
  sequences never cut.

## vsearch

redig is `--cut` with `--cut_pattern` as `--pattern`, `--fastaout_rev` as
`--rev-out`, `--fastaout_discarded` as `--uncut` and
`--fastaout_discarded_rev` as `--rev-uncut`. Its outputs are unwrapped.
vsearch's are `conv --fasta` of them.

## Examples

```console
$ cat digest.fa
>d1
AAACCCGAATTCTTGCAG
>d2
ACGTACGTACGT
$ redig --pattern G^AATT_C --rev-out rev.fa --uncut uncut.fa -i digest.fa
>d1
AAACCCG
>d1
AATTCTTGCAG
$ cat rev.fa
>d1
AATTCGGGTTT
>d1
CTGCAAG
$ cat uncut.fa
>d2
ACGTACGTACGT
```

## Exit status

0 on success, 1 for invalid input (including FASTQ), 2 for invalid
arguments, such as a pattern without both cut sites.

```console
$ redig --pattern GAATTC -i digest.fa; echo "exit $?"
error: invalid value 'GAATTC' for '--pattern <PATTERN>': no forward cut site (^) in the pattern

For more information, try '--help'.
exit 2
```

## See also

[subseq](subseq.md), [conv](conv.md)
