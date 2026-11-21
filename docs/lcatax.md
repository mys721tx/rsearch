# lcatax

Write the last common ancestor of each query's top hits.

## Synopsis

```
lcatax HITS --db FILE [-o FILE] [--cutoff SHARE]
```

## Description

lcatax reads UC hits and the database they were found in (`--db`: FASTA,
FASTQ or UDB). For each query of the hits it writes its header, a tab, and
the `tax=` levels its targets agree on. It goes level by level, from the
first, while at least `--cutoff` of the hits agree on a name, as vsearch's
majority vote picks it. A query without hits gets an empty taxonomy.

The hits are those to be counted. Limit them in the search
(`--maxhits`, `--top-hits-only`), as vsearch's `--maxhits` and
`--top_hits_only` limit the hits its `--lcaout` reads.

## Options

- `HITS`: UC file of hits. Required.
- `--db FILE`: Database the hits were found in: FASTA, FASTQ or UDB file.
  Required.
- `-o`, `--out FILE`: Output file. Default: standard output.
- `--cutoff SHARE`: Share of the hits that must agree on a level, above 0.5
  and at most 1. Default: 1.

## vsearch

lcatax is the `--lcaout` output of vsearch's search commands, with
`--lca_cutoff` as `--cutoff`.

## Examples

```console
$ dbsrch -i queries.fa --db db.fa --id 0.85 --maxaccepts 0 > hits.uc
$ lcatax hits.uc --db db.fa
q1;sample=S1;size=3	d:Bacteria,p:Firmicutes,g:Bacillus
q2;sample=S1;size=1	d:Bacteria,p:Proteobacteria,g:Escherichia,s:Escherichia_coli
q3;sample=S2;size=2	d:Bacteria,p:Firmicutes,g:Bacillus
q4;sample=S2;size=1	
q5;sample=S2;size=1	
```

## Exit status

0 on success, 1 for invalid input, 2 for invalid arguments, such as a
`--cutoff` of 0.5 or less.

## See also

[dbsrch](dbsrch.md), [sintax](sintax.md), [otutab](otutab.md)
