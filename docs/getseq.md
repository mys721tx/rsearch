# getseq

Select records by their headers or by hits.

## Synopsis

```
getseq [-i FILE] [-o FILE] (--label LABEL | --labels FILE | --word WORD |
       --words FILE | --hits FILE [--targets]) [--substr] [--field NAME]
       [--invert] [--discarded FILE]
```

## Description

getseq writes the records it selects, unchanged, in the input format, FASTA
unwrapped; with `--invert`, the others. `--discarded` takes the records not
written. One way of selecting is given:

- A label (`--label`, or a file of labels one per line with `--labels`)
  matches the whole header, ignoring case, or any part of it with
  `--substr`.
- A word (`--word`, or a file of words with `--words`) matches in the header
  between characters that are not letters or digits, case included; with
  `--field NAME`, it matches the value of `NAME=` between semicolons.
- `--hits` takes a UC file of hits, from [dbsrch](dbsrch.md),
  [exsrch](exsrch.md), [allpr](allpr.md) or the clustering tools. The input
  is the queries of the search, matched in order to the queries of the UC
  file by label. A query is selected when it has a hit, and not when it has
  an N record. With `--targets`, the input is the database instead. A
  record is then selected when some hit names it, by its number.

Each line of a labels or words file is read whole.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the selected records. Default:
  standard output.
- `--label LABEL`: Select the records with this header.
- `--labels FILE`: Select the records with a header listed in FILE.
- `--word WORD`: Select the records whose header holds this word.
- `--words FILE`: Select the records whose header holds a word listed in
  FILE.
- `--hits FILE`: Select the queries with hits in this UC file, matched in
  order.
- `--targets`: With `--hits`, select the database records hit, by their
  number.
- `--substr`: Match labels as any part of the header. Not with `--word` or
  `--words`.
- `--field NAME`: Match words as the value of `NAME=` between semicolons.
  Not with `--label` or `--labels`.
- `--invert`: Write the records that are not selected instead.
- `--discarded FILE`: Output file of the records not written to the output.

## vsearch

getseq stands for `--fastx_getseq` (`--label`, `--label_substr_match`),
`--fastx_getseqs` (`--labels`, `--label_word`, `--label_words`,
`--label_field`), the selection half of `--fastx_getsubseq`, and the
`--matched`, `--notmatched`, `--dbmatched` and `--dbnotmatched` outputs of
the search commands. `--hits` is `--matched`, `--hits --invert`
`--notmatched`, `--hits --targets` `--dbmatched` and `--hits --targets
--invert` `--dbnotmatched`. vsearch's `--sizeout` on `--dbmatched` is
[relab](relab.md). vsearch counts every accepted hit for `--dbmatched`. The
search must therefore write them all (no `--maxhits` or `--top-hits-only`).
vsearch reads a labels or words file in pieces of 1023 bytes, and a longer
line becomes several labels. getseq reads each line whole.

## Examples

```console
$ getseq --label S3\;SIZE=2 -i seqs.fa
>s3;size=2
TTGACCGTAGCA
$ getseq --labels labels.txt --invert -i seqs.fa
>s1;size=3
ACGTTGCAAGGT
>s2
acgttgcaaggt
>s4
ACCTTGCAAGGTAAC
$ getseq --word 2 --field size -i seqs.fa
>s3;size=2
TTGACCGTAGCA
```

The queries with hits, and the database records hit, as vsearch's
`--usearch_global queries.fa --db db.fa --id 0.9 --matched m.fa --dbmatched
d.fa`:

```console
$ dbsrch -i queries.fa --db db.fa --id 0.9 > hits.uc
$ getseq --hits hits.uc -i queries.fa | grep '>'
>q1;sample=S1;size=3
>q2;sample=S1;size=1
>q3;sample=S2;size=2
$ getseq --hits hits.uc --targets -i db.fa | grep '>'
>ref1;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis
>ref2;tax=d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_cereus
>ref3;tax=d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli
```

## Exit status

0 on success, 1 for invalid input (including a query of the hits missing
from the input), 2 for invalid arguments.

## See also

[subseq](subseq.md), [topn](topn.md), [dbsrch](dbsrch.md),
[exsrch](exsrch.md), [clspl](clspl.md)
