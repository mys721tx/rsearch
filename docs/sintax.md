# sintax

Classify sequences by taxonomy with SINTAX.

## Synopsis

```
sintax --db FILE [-i FILE] [-o FILE] [--word-len N] [--strand plus|both]
       [--cutoff SHARE] [--random] [--seed SEED] [--report FILE]
```

## Description

sintax classifies each query as vsearch `--sintax` does. The database
(`--db`: FASTA, FASTQ or UDB, with `;tax=` annotations) is read whole and its
words of `--word-len` indexed. The queries are then read in turn. A hundred
times, 32 words drawn at random from the query pick the database sequence
sharing the most of them (ties go to the shortest, then first, sequence, or
at random with `--random`). With fifty picks or more, each level of the
`tax=` annotations takes the name most picks agree on.

The table has a line per query: its header, the names with the share of
picks agreeing (as `d:Bacteria(1.00),...`), the strand, and with `--cutoff`
the names whose share reaches it. A query not classified has empty fields.
Lower-case residues are masked, in the queries and the database.

With a seed, the result is the one vsearch gives with `--randseed`. Without
one, or with 0, the seed comes from the system, and the result can change
from run to run.

## Options

- `-i`, `--in FILE`: Input FASTA or FASTQ file of queries. Default:
  standard input.
- `-o`, `--out FILE`: Output table. Default: standard output.
- `--db FILE`: Database: FASTA, FASTQ or UDB file with `;tax=` annotations.
  Required.
- `--word-len N`: Length of the words, 3 to 15; a UDB database has its own.
  Default: 8.
- `--strand plus|both`: Search the plus strand, or both. Default: plus.
- `--cutoff SHARE`: Also write the names agreed on by at least this share
  of the picks.
- `--random`: Break ties at random, rather than for the shortest, then
  first, sequence.
- `--seed SEED`: Seed of the random numbers. Default: one from the system.
- `--report FILE`: Output file of the number of queries classified.

## vsearch

sintax is `--sintax` with `--tabbedout` as `-o`, `--randseed` as `--seed`,
`--sintax_cutoff` as `--cutoff`, `--sintax_random` as `--random` and
`--wordlength` as `--word-len`. The count of queries classified, which
vsearch prints, goes to `--report`. vsearch masks the lower case of its
database (for `--dbmask dust` and `soft` alike, without running DUST) and
never masks its queries. Its result is therefore sintax with `lcmask
--clear` on the queries. `--dbmask none` is `lcmask --clear` on the
database too. Its
database length limits are `sqfilt --min-len 32 --max-len 50000` on the
database.

## Examples

```console
$ sintax -i queries.fa --db db.fa --seed 1 --cutoff 0.8 --strand both | cut -f 1,3,4 | head -4
q1;sample=S1;size=3	+	d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_subtilis
q2;sample=S1;size=1	+	d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli
q3;sample=S2;size=2	+	d:Bacteria,p:Firmicutes,g:Bacillus,s:Bacillus_cereus
q4;sample=S2;size=1	-	d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli
```

vsearch's `--sintax queries.fa --db db.fa --randseed 1 --tabbedout t.tsv`:

```console
$ sqfilt --min-len 32 --max-len 50000 -i db.fa > dbf.fa
$ lcmask --clear -i queries.fa | sintax --db dbf.fa --seed 1 --report report.txt | cut -f 1,3
q1;sample=S1;size=3	+
q2;sample=S1;size=1	+
q3;sample=S2;size=2	+
q4;sample=S2;size=1	
q5;sample=S2;size=1	
$ cat report.txt
Classified 3 of 5 sequences (60.00%)
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments.

## See also

[lcatax](lcatax.md), [dbsrch](dbsrch.md), [mkudb](mkudb.md),
[lcmask](lcmask.md)
