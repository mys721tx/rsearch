# conv

Change the format of records and the encoding of their quality scores.

## Synopsis

```
conv [-i FILE] [-o FILE] [--ascii 33|64] [--qmin Q] [--qmax Q] [--solexa]
     [--ascii-out 33|64] [--qmin-out Q] [--qmax-out Q] [--clip]
     [--fasta [--fasta-width N]]
```

## Description

conv reads FASTQ, FASTA, SFF or UDB, recognised from the start of the
input, and writes FASTQ, or FASTA with `--fasta`. It is the only tool that
changes the format of records or wraps FASTA lines. Every other tool writes
records in the format it read, FASTA unwrapped. `... | conv --fasta` thus
turns any pipeline's output into the FASTA a vsearch command writes by
default.

- FASTQ scores are read with the input offset (`--ascii`) and checked
  against `--qmin` and `--qmax`; with `--solexa` they are converted from the
  Solexa scale to Phred. They are then limited to `--qmin-out` and
  `--qmax-out` and written with the output offset (`--ascii-out`).
- FASTA sequences get the highest output score at every base (41 unless
  `--qmax-out` says otherwise).
- SFF reads keep their Phred scores; the bases outside the clip points are
  written in lower case, or removed with `--clip`.
- UDB files are read as FASTA: their sequences and headers as stored.
- With `--fasta`, the scores of any input are checked, then dropped.
  Sequence lines are wrapped at `--fasta-width` (80, vsearch's default; 0
  for one line).

Headers are written unchanged.

## Options

- `-i`, `--in FILE`: Input FASTQ, FASTA, SFF or UDB file. Default: standard
  input.
- `-o`, `--out FILE`: Output FASTQ file, or FASTA with `--fasta`. Default:
  standard output.
- `--ascii 33|64`: ASCII offset of the input quality scores. Default: 33.
  FASTQ input only.
- `--qmin Q`: Lowest input score accepted. Default: 0, or -5 with
  `--solexa`. FASTQ input only.
- `--qmax Q`: Highest input score accepted. Default: 126 minus `--ascii`.
  FASTQ input only.
- `--solexa`: The input scores are on the Solexa scale; convert them to
  Phred. Needs `--ascii 64`.
- `--ascii-out 33|64`: ASCII offset of the output quality scores.
  Default: 33.
- `--qmin-out Q`: Lowest score written; lower ones are raised to it.
  Default: 0.
- `--qmax-out Q`: Highest score written; higher ones are lowered to it.
  Default: 126 minus `--ascii-out`, or 41 for FASTA input.
- `--clip`: Remove the bases outside the clip points of SFF reads, instead
  of writing them in lower case.
- `--fasta`: Write FASTA instead of FASTQ.
- `--fasta-width N`: Wrap FASTA sequence lines at N residues, 0 for no
  wrapping. Default: 80. Needs `--fasta`.

An option that needs quality scores (`--ascii`, `--qmin`, `--qmax`,
`--solexa`, and `--qmin-out` for FASTA input) is refused for input that has
none. The output score options are refused with `--fasta`.

## vsearch

conv stands for `--fastq_convert`, `--fasta2fastq`, `--sff_convert` and
`--udb2fasta`, and for the `--fastaout` and `--fasta_width` of every vsearch
command. Its options are vsearch's `--fastq_ascii`, `--fastq_qmin`,
`--fastq_qmax`, `--fastq_asciiout`, `--fastq_qminout`, `--fastq_qmaxout` and
`--sff_clip`; `--solexa` replaces vsearch's choice of Solexa input by
`--fastq_ascii 64` with `--fastq_qmin -5`. Unlike vsearch, conv never cuts a
header at its first space or tab, as vsearch does unless given
`--notrunclabels`.

## Examples

FASTQ to FASTA, wrapped at 10 residues:

```console
$ conv --fasta --fasta-width 10 -i reads.fq
>r1
ACGTACGTAG
CTAGCTAGGA
TCGATCGA
>r2
TTGACCGTAG
CAAGCTTAGC
NAGCTAAT
>r3
GGCATCGATC
GA
>r4
ACGTACGTAG
CTAGCTAGGA
TCGATCGA
```

FASTA to FASTQ, and FASTQ scores to offset 64, capped at 30:

```console
$ conv -i seqs.fa | head -4
@s1;size=3
ACGTTGCAAGGT
+
JJJJJJJJJJJJ
$ conv --ascii-out 64 --qmax-out 30 -i reads.fq | head -4
@r1
ACGTACGTAGCTAGCTAGGATCGATCGA
+
^^^^^^^^^^^]\[YXWVUTSRQPONML
```

A UDB file back to FASTA:

```console
$ mkudb -i seqs.fa -o seqs.udb
Warning: 1 of 5 sequences give no word to index: they are shorter than the word length (8), or each of their words holds a lower-case residue or a symbol other than A, C, G, T or U
$ conv --fasta -i seqs.udb
>s1;size=3
ACGTTGCAAGGT
>s2
acgttgcaaggt
>s3;size=2
TTGACCGTAGCA
>s4
ACCTTGCAAGGTAAC
>s5
GGGAAATTTCCC
```

## Exit status

0 on success, 1 for invalid input (a malformed record, a score out of
range, an unreadable SFF or UDB file), 2 for invalid arguments, such as an
option that needs scores on FASTA input:

```console
$ conv --ascii 64 -i seqs.fa; echo "exit $?"
error: --ascii cannot be used with FASTA input

Usage: conv [OPTIONS]

For more information, try '--help'.
exit 2
```

## See also

[sqtrim](sqtrim.md), [sqfilt](sqfilt.md), [relab](relab.md),
[mkudb](mkudb.md)
