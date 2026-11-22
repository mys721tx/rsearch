# dbornt

Orient sequences by database words; write them forward.

## Synopsis

```
dbornt --db FILE [-i FILE] [-o FILE] [--word-len N] [--not-oriented FILE]
       [--tabbedout FILE] [--report FILE] [-t N]
```

## Description

dbornt finds the strand of each query from the words it shares with a
database of sequences in a known orientation, as vsearch `--orient` does.
It first checks the queries' format and creates its outputs. It then reads
the whole database (`--db`: FASTA, FASTQ or UDB) and indexes its words of
`--word-len`. Last, it reads the queries in turn. Each distinct word of a
query votes for a strand when more than eight times as many database
sequences hold it on that strand as on the other. A query is oriented when
one strand has a vote and four times the votes of the other.

Queries on the plus strand are written as they are, those on the minus
strand reverse-complemented. Both go to `-o` in their format, FASTA
unwrapped. The others go to `--not-oriented`. `--tabbedout` writes a line
per query: its header, the strand (`+`, `-` or `?`) and the votes for each
strand. `--report` writes the counts. Lower-case residues are masked, in the
queries and the database.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file. Default: standard input.
- `-o`, `--out FILE`: Output file of the oriented sequences. Default:
  standard output.
- `--db FILE`: Database: FASTA, FASTQ or UDB file. Required.
- `--word-len N`: Length of the words, 3 to 15; a UDB database's own length
  wins, with a warning when they differ. Default: 12.
- `--not-oriented FILE`: Output file of the sequences not oriented.
- `--tabbedout FILE`: Output table: header, strand (`+`, `-` or `?`) and
  the votes for each strand.
- `--report FILE`: Output file of the summary.
- `-t`, `--threads N`: Orient on N threads, 0 for one per logical core.
  The queries go to the threads 1,000 at a time, and the outputs are the
  same, in the same order, on any number. Default: 0.

## vsearch

dbornt is `--orient`, with `--notmatched` as `--not-oriented`,
`--wordlength` as `--word-len`, and the summary vsearch prints written to
`--report`. Its `--fastaout` from FASTQ input is `conv --fasta` after it.
vsearch reads its database as it is, unlike other commands. `--dbmask soft`
is the database without [lcmask](lcmask.md), `--dbmask none` `lcmask
--clear`, and the default `--dbmask dust` `lcmask`. It masks the lower case
of its queries without running DUST. dbornt therefore takes them as they
are. `--qmask none` is `lcmask --clear` on them. Its database length limit
is `sqfilt --max-len 50000`.

vsearch orients on one thread. dbornt orients on several and writes the
same as on one.

## Examples

```console
$ dbornt -i queries.fa --db db.fa --tabbedout strands.tsv --report report.txt | grep '>'
>q1;sample=S1;size=3
>q2;sample=S1;size=1
>q3;sample=S2;size=2
>q4;sample=S2;size=1
$ cat strands.tsv
q1;sample=S1;size=3	+	137	0
q2;sample=S1;size=1	+	149	0
q3;sample=S2;size=2	+	149	0
q4;sample=S2;size=1	-	0	149
q5;sample=S2;size=1	?	0	0
$ cat report.txt
Forward oriented sequences: 3 (60.00%)
Reverse oriented sequences: 1 (20.00%)
All oriented sequences:     4 (80.00%)
Not oriented sequences:     1 (20.00%)
Total number of sequences:  5
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[revcomp](revcomp.md), [lcmask](lcmask.md), [mkudb](mkudb.md),
[conv](conv.md)
