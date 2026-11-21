# derep

Group identical sequences and count them.

## Synopsis

```
derep [-i FILE] [-o FILE] [-s | --strand] [--same-label] [--prefix]
      [--unsorted] [--uc FILE] [--tabbedout FILE] [--no-sizein]
      [--no-sizeout] [--ascii 33|64] [--qmin Q] [--qmax Q]
      [--ascii-out 33|64] [--qmin-out Q] [--qmax-out Q] [--best-qual]
```

## Description

derep groups the sequences that are identical, ignoring case and taking U
for T. It writes each group once, as its first sequence and header with the
group's abundance appended as `;size=`. The abundance is the sum of the
members' `;size=`, 1 for a member without one (each member counts 1 with
`--no-sizein`).

- With `--strand`, a sequence also joins the group of its reverse
  complement.
- With `--same-label`, only sequences with identical headers are grouped.
- With `--prefix` (FASTA only), each sequence joins the longest sequence it
  is a prefix of.
- FASTQ scores are merged: each position gets the average of the members'
  error probabilities, or with `--best-qual` the best score.

Groups are written by decreasing abundance, then header, then first
appearance; with `--unsorted`, in order of first appearance. Records are
written in the input format, FASTA unwrapped. `--uc` writes the groups in
UC format (an S record for each group, an H record for each further member,
a C record for each group). `--tabbedout` writes one line per member: its
header, the group's header, the group's number, the member's number in the
group, the group's member count and the group's header again.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file, in the input format. Default: standard
  output.
- `-s`, `--strand`: Also group reverse-complemented duplicates.
- `--same-label`: Group only duplicates that also have identical headers.
- `--prefix`: Group sequences into the longest sequence they are a prefix
  of. FASTA only; not with `--strand` or `--same-label`.
- `--unsorted`: Write the groups in order of first appearance instead of
  decreasing abundance.
- `--uc FILE`: Output file of the groups in UC format.
- `--tabbedout FILE`: Output file of every group's members, tab-separated.
- `--no-sizein`: Count each sequence once, ignoring its `;size=`.
- `--no-sizeout`: Do not append the abundance to the headers written.
- `--ascii 33|64`: ASCII offset of the input quality scores. Default: 33.
- `--qmin Q`: Lowest input score accepted. Default: 0.
- `--qmax Q`: Highest input score accepted. Default: 126 minus `--ascii`.
- `--ascii-out 33|64`: ASCII offset of the output quality scores.
  Default: 33.
- `--qmin-out Q`: Lowest score written. Default: 0.
- `--qmax-out Q`: Highest score written. Default: 126 minus `--ascii-out`.
- `--best-qual`: Keep the best score of the merged sequences instead of
  averaging their error probabilities.

The score options are refused for FASTA input.

## vsearch

derep stands for `--fastx_uniques`, `--derep_fulllength`, `--derep_id`
(`--same-label`), `--derep_prefix` (`--prefix`) and `--derep_smallmem`
(`--unsorted`); `--fastq_qout_max` is `--best-qual`. It differs in what it
does by default:

- It reads `;size=` and writes it unless told not to, where vsearch does
  neither unless given `--sizein` and `--sizeout`. vsearch's
  `--derep_fulllength in --output out` is `derep --no-sizein --no-sizeout`.
- It never relabels and never filters. `--relabel` and its kin are
  [relab](relab.md) after it. `--minuniquesize`, `--maxuniquesize` and
  `--minseqlength` are [sqfilt](sqfilt.md), and `--topn` is
  [topn](topn.md).
  As it never relabels, its `--tabbedout` always shows the first member's
  header, which differs from vsearch's only when vsearch relabels.
- Headers are whole, as with vsearch's `--notrunclabels`. Abundances,
  labels compared by `--same-label`, `--uc` and `--tabbedout` all use the
  whole header. A `;size=` followed by a description is therefore not an
  abundance.

## Examples

```console
$ cat seqs.fa
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
$ derep -i seqs.fa
>s1;size=4
ACGTTGCAAGGT
>s3;size=2
TTGACCGTAGCA
>s4;size=1
ACCTTGCAAGGTAAC
>s5;size=1
GGGAAATTTCCC
$ derep --unsorted --no-sizein -i seqs.fa --uc groups.uc
>s1;size=2
ACGTTGCAAGGT
>s3;size=1
TTGACCGTAGCA
>s4;size=1
ACCTTGCAAGGTAAC
>s5;size=1
GGGAAATTTCCC
$ cat groups.uc
S	0	12	*	*	*	*	*	s1;size=3	*
H	0	12	100.0	+	0	0	*	s2	s1;size=3
S	1	12	*	*	*	*	*	s3;size=2	*
S	2	15	*	*	*	*	*	s4	*
S	3	12	*	*	*	*	*	s5	*
C	0	2	*	*	*	*	*	s1;size=3	*
C	1	1	*	*	*	*	*	s3;size=2	*
C	2	1	*	*	*	*	*	s4	*
C	3	1	*	*	*	*	*	s5	*
```

vsearch's `--fastx_uniques seqs.fa --minuniquesize 2 --relabel U --sizeout
--fastaout out.fa`:

```console
$ derep --no-sizein -i seqs.fa | sqfilt --min-size 2 | relab --relabel U --sizeout | conv --fasta
>U1;size=2
ACGTTGCAAGGT
```

## Exit status

0 on success, 1 for invalid input (including an invalid `;size=`), 2 for invalid
arguments, such as `--prefix` on FASTQ input.

```console
$ derep --prefix -i reads.fq; echo "exit $?"
error: --prefix cannot be used with FASTQ input

Usage: derep [OPTIONS]

For more information, try '--help'.
exit 2
```

## See also

[rerep](rerep.md), [sqfilt](sqfilt.md), [relab](relab.md),
[sortby](sortby.md), [clust](clust.md)
