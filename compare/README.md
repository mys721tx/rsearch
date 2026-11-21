# Comparing rsearch with vsearch

Untracked helpers for checking that rsearch's tools give vsearch's
results. They form the Rust crate `rscompare`. It is a member of the
workspace. It translates a vsearch command line into the rsearch
pipeline that should give the same results and runs both. It then
compares their outputs, exit codes, error messages, and warnings. Its
dependencies are `clap` and `serde_json`. It calls `gzip`, `bzip2`, and,
for `fetch`, `curl` and `md5`/`md5sum`.

## Build

From the repository root, build rsearch's tools and rscompare into
`target/release`, then vsearch from `../vsearch` into
`compare/vsearch-build/vsearch`:

```sh
cargo build --release
./target/release/rscompare build-vsearch
```

`build-vsearch [SOURCE] [OUT]` compiles vsearch's sources directly (no
autotools). It chooses the SIMD backends of `src/Makefile.am` for the
machine's architecture and the system libraries for its OS. `--arch`,
`--os`, or `VSEARCH_ARCH`/`VSEARCH_OS` override them, and `--dry-run`
prints the choice. It defines `HAVE_ZLIB_H` and `HAVE_BZLIB_H`, as
configure would, so that vsearch reads gzip and bzip2 files. It opens
`libz` and `libbz2` at run time. Without them it fails with "Files
compressed with gzip are not supported". `CXX` and `CC` are honoured.

## Compare

From the repository root, `list` shows the sources, their default seeds,
and what they cover. `run` runs sources at their default seeds, or at
those of `--seeds`, and `summary` shows the outcomes, the coverage, and
the caveats of a results file. `translate` prints the pipeline of a
vsearch command line, and `cargo test -p rscompare` runs the unit and
translator tests. `cargo run -p rscompare --release --` can stand for
`./target/release/rscompare`:

```sh
./target/release/rscompare list
./target/release/rscompare run derep clust
./target/release/rscompare run all --out compare/results/run.jsonl
./target/release/rscompare run dbsrch --seeds 1-50,77 -j 8
./target/release/rscompare summary compare/results/run.jsonl
./target/release/rscompare translate -- --cluster_fast in.fa --id 0.97 --centroids c.fa
cargo test -p rscompare
```

`run` prints a `DIFF` line for each failing case (the first 200). It
then prints the tally per source of `pass`, `fail`, `skip`, and `error`
(a case the translator refused, or a harness failure). It exits 1 if any
case failed.

`--vsearch` and `--bin` name the vsearch binary and the rsearch build
(defaults `compare/vsearch-build/vsearch` and `target/release`).
`--jobs` sets the cases run at once (threads; each tool stays
single-threaded). `--tmp` sets the working directory
(`$TMPDIR/rscompare`). `--keep` keeps the files of passing cases too.
The files of a failing case are always kept in
`$TMPDIR/rscompare/<source>-<seed>/c<index>/`. vsearch's outputs are in
`v/` and rsearch's in `r/`.

`all` runs every source but two. The UDB files of `mkudb-long` reach 4
GB, and `datasets` needs the downloads (below). `mkudb-long` runs on its
own, one case at a time:

```sh
./target/release/rscompare run mkudb-long -j 1
```

### Results

`--out` appends one JSON object a case:

| key | what |
|---|---|
| `source`, `seed`, `index`, `case` | where the case comes from |
| `kind` | `translated`, or `probe` for a case that gives its rsearch steps itself |
| `vsearch` | vsearch's argv |
| `rsearch` | the rsearch pipeline, a line a chain of steps |
| `outcome` | `pass`, `fail`, `skip`, or `error`. The last two have a `reason` |
| `diffs` | what differed: exit code, error, a warning, or the first different line of an output |
| `exit`, `stderr` | both programs' exit codes, and their error messages on a failure |
| `failed_step` | the first rsearch step that failed |
| `coverage` | the case's parameters, vsearch's exit code and error, and counts from vsearch's outputs (UC record types, records, table rows) |
| `caveats` | what the translation does not reproduce (below) |
| `seconds`, `wall` | time in vsearch, in the rsearch steps, and the whole case |

`summary` tallies the outcomes, the coverage values, and the caveats of
one or more result files so that two runs can be compared.

## The translator

`src/translate.rs` can be used alone:

```rust
let p = rscompare::translate::translate(&["--fastx_uniques", "in.fq", "--sizeout", "--fastqout", "o.fq"], "work")?;
p.steps      // [Step { tool, args, stdin, stdout, system }] in the order they run
p.outputs    // [Output { option, vsearch, rsearch, comparator }]
p.caveats    // what the pipeline does not reproduce
println!("{}", p.render());
```

It holds each rule of the translation once. The rules come from the
pipelines of `README.md`, the tool pages in `docs/`, and the option map
of `docs/vsearch.md`. They are:

- option names;
- `sqtrim` then `sqfilt` for the filters;
- `sqfilt` for vsearch's length filters (`--minseqlength`, 32 or 1 by
  command, `--maxseqlength 50000`), `--minsize`, `--minuniquesize`, and
  the unmasked shares;
- `lcmask` for `--qmask`, `--dbmask`, and `--hardmask` by each command's
  rules (searches, `--makeudb_usearch`, `--orient`, `--sintax`);
- `sortby` before `--cluster_fast`, `--cluster_size`,
  `--cluster_unoise`, and de novo chimeras;
- `--maxrejects 8` for `--cluster_fast`;
- `--sizein`/`--sizeout` against rsearch's defaults;
- `topn` and `relab` for the label options, and `conv --fasta` at
  `--fasta_width`;
- `hitfmt`, `getseq --hits`, `otutab`, `lcatax`, `clspl`, and `clmsa` on
  the hits;
- a search without `--maxhits` for `--dbmatched`, one with `--maxhits 1`
  for `--uc` without `--uc_allhits`;
- side outputs (`clust --hits`, `dbornt --not-oriented`, `redig`'s extra
  files) only when vsearch asks for them;
- which rsearch output stands for each vsearch output. `--log` is the
  report of `--fastq_chars`, `--fastq_stats`, `--fastq_mergepairs`,
  `--orient`, `--sintax`, `--udbinfo`, and `--udbstats`, and otherwise a
  log that is not compared.

`execute::run_steps` runs the steps one after the other through files.
It is a shell pipeline without the concurrency. As in a pipeline, a
failed step's output still streams into the next step on its standard
input. A step that takes it as a file argument (a database, a hit file)
does not run. The exit code and error are the first failing step's. For
the benchmarks, `execute::run_piped` runs the same steps as a real
pipeline.

Compressed input is decompressed first. vsearch knows a gzip or bzip2
file by its magic number, and so does the translator. The main input,
`--reverse`, and `--db` go through `gzip -dc` or `bzip2 -dc` into the
pipeline. `--gzip_decompress` and `--bzip2_decompress` matter to vsearch
only for a pipe, and it reads an uncompressed file as it is with them.
The translator uses `gzip -dc` because macOS's `zcat` expects `.Z`
files. A failing decompressor counts as exit code 1. That is vsearch's
code for unreadable input. What it writes before failing is its own.
From a file cut short, macOS `gzip -dc` writes nothing. vsearch writes
the records it read.

An option it cannot translate is a `TranslationError` naming it (for
example `--alnout` with clustering, `--fastq_minlen 0`, compressed
standard input, and `--clusters` at a width other than 0). A default
behaviour of vsearch that the pipeline does not reproduce is a caveat.
Each case records its caveats:

- `--minseqlength 0` keeps empty sequences. `sqfilt` cannot keep them,
  and `--maxseqlength` is then not applied.
- Without `--notrunclabels`, vsearch cuts labels at spaces. rsearch
  keeps them whole.
- Without `--threads 1`, the multithreaded commands may order their
  outputs otherwise.
- Without `--randseed`, random draws differ.
- `sqfilt` checks `;size=`. Dereplication reads it only with `--sizein`.
- `getseq --targets` counts weak hits, and vsearch's `--dbmatched`
  counts accepted hits only. `clspl` writes the abundances of
  `--clusters` with `--sizeout` as they are. vsearch's `--search_exact`
  writes masked sequences.

## Benchmarks

`benches/<tool>/main.rs` (`cargo bench -p rscompare`) times vsearch against
rsearch. Each benchmark times vsearch commands (`vsearch/<case>`). For each,
it also times the pipeline the translator gives for that command
(`pipeline/<case>`: `rscompare::translate` then `execute::run_piped`). What
is timed against vsearch is thus by construction what the comparison checks.
Where a tool's output can be compared directly, the benchmark first checks
that rsearch and vsearch write the same.

The benchmarks run on the input files of rsearch's own benchmark of the
same tool. Both include its generator (`benches/<tool>/input.rs` at the
root). They use the same criterion group names, and `target/criterion`
holds the rsearch, pipeline, and vsearch timings of a case side by side.
rsearch's benchmarks time rsearch alone.

vsearch is `VSEARCH`, or else the build of `rscompare build-vsearch`.
Without either, the benchmarks time nothing. The pipelines run rsearch's
release tools. Build them first:

```sh
cargo build --release
cargo bench -p rscompare --bench derep
```

The vsearch command lines keep the options that leave vsearch's default
filters, masking and wrapping out where they change nothing on the
benchmark's input (`--minseqlength 0`, `--qmask soft`, and
`--fasta_width 0`). They run vsearch with `--threads 1` so that the
rsearch timings stay comparable.

## Comparators

`src/comparators.rs` holds the normalisations. It names each once:

| comparator | for |
|---|---|
| `bytes` | everything else |
| `uc-fast` | `--cluster_fast`'s `=` for an alignment perfect but for terminal gaps. clust writes the CIGAR in full there |
| `biom` | BIOM without its id, program, and date lines. Those lines name otutab's output, rsearch, and the time |
| `log-report` | a report in vsearch's `--log` without the log's own lines |
| `orient-log`, `sintax-log`, `udbinfo-log`, `udbstats-log` | the parts of `--log` those tools write to `--report` or `-o` |
| `alnout`, `sam` | `--alnout` without its opening lines, SAM without `@PG` and the database's file name |
| `cluster-files` | the files of `--clusters` by number |
| `upper-seq`, `fastq-to-fasta` | vsearch's output upper-cased or in FASTA for two probes |

Errors are compared by kind: `exit` (exit code only), `exact` (vsearch's
`Fatal error: X` against rsearch's `Error: X`), or a mapping of messages
that differ in wording (`abundance`, `joinpa`, `mergpa`, `mkudb`,
`mkudb-exact`, `udb`, `invalid`). The warnings compared are `stripped`,
`missing-abundance`, `nowords`, and `wordlength-adjusted`. rsearch's
exit code 2 for invalid arguments stands for vsearch's 1 where a case
says so.

## Sources and generators

`src/sources/` makes the cases. `src/generators/` writes their inputs
from a seed (`src/rng.rs`, SplitMix64). Every case is reproducible from
its seed on any platform.

- `preprocess.rs` and `analysis.rs`: a source for each tool, or for a
  use of one (`getseq-hits`, `conv-udb`). Each has its own option sets
  and seeds. Where a source applies an rsearch tool to vsearch's own
  hits (hitfmt, otutab, lcatax, getseq `--hits`, and clspl), the case is
  vsearch's command whole. The translated pipeline then runs the search
  too.
- `fuzz.rs`: `fuzz` (malformed FASTA and FASTQ against about 40 vsearch
  commands) and `fuzz-report` (the quality reports). Where the
  translated pipeline starts with another tool than the one a case
  fuzzes, a probe runs that tool's own command too.
- `compressed.rs`: `compressed` gzips or bzip2s the generators' files at
  random (seeds 1 to 100). It covers streaming and whole-file tools, a
  FASTA output at 80, pairs, and a search with a compressed database. It
  also covers fuzzed records inside compressed files and files cut short
  (exit code only).
- `datasets.rs`: the real datasets of `datasets.json` (below). Seed N is
  the manifest's Nth entry. An entry's command lines are read from
  `"vsearch"`, `"commands"`, or `"vsearch_commands"` (strings or argv
  lists, without the program). `{dir}` in them stands for the dataset's
  directory, and relative input names are taken from it. Output options
  are moved to the case's output directory. `--threads 1`, `--quiet`,
  and `--notrunclabels` are added. An entry whose files are missing is
  skipped.

A source is a function from a directory, a seed, and the environment to
`Case`s (`src/cases.rs`). A case holds vsearch's argv with outputs under
`@OUT@` and the error kind and warnings to compare. It also says whether
outputs are compared when vsearch fails. It may hold a probe or a skip
rule.

The 37 generators write each source's inputs from its seed: the edge
cases, the option draws, the files, and the `params`.

## Real datasets

```sh
./target/release/rscompare fetch --list
./target/release/rscompare fetch --max-reads 5000 --vsearch compare/vsearch-build/vsearch
./target/release/rscompare run datasets
```

`fetch [ID ...]` downloads the datasets of `datasets.json` (and the
references they use) into `$RSCOMPARE_DATA/<id>/` (default
`$TMPDIR/rsearch-data/<id>/`). It resumes partial downloads and checks
sizes and MD5s. It then decompresses the files and cuts each FASTQ to
its first N reads with `--max-reads`. Last, it runs the entry's
`prepare` steps (`cat`, `fasta`, `derep`, and `makeudb` with
`--vsearch`). The working files are rebuilt only when the downloads,
`--max-reads`, or the steps change. See `docs/datasets.md`.

## Time two builds

```sh
./target/release/rscompare time BASELINE_DIR target/release
```

Runs the two builds alternately on the benchmark inputs in `target/tmp`
(written by `cargo bench -- --test`) and prints the median time of each.
`--runs` sets the runs of each build (15 by default), and `--inputs` the
inputs. A baseline from another commit, such as `HEAD`, builds from its
files into `/tmp/base-target/release`. That directory is the
`BASELINE_DIR` above:

```sh
mkdir -p /tmp/base
git archive HEAD | tar -x -C /tmp/base
(cd /tmp/base && CARGO_TARGET_DIR=/tmp/base-target cargo build --release)
```

## Runs

- First full run (2026-11-20, rsearch `26fdb8c`, vsearch `e1551256`;
  `results/final-2026-11-20.jsonl` and its `.summary.txt`). It ran every
  source at its default seeds. 71,570 cases passed and none failed. The
  sources' own rules skipped 18 (14 `conv-udb` files whose magic number
  is spoiled, and 4 `udbstat` files on which vsearch aborts).
  `compressed` was added later and passed 1,200 cases
  (`results/compressed-2026-11-20.jsonl`).
- A run on the present generators (`results/rust-gen.jsonl`, 44 sources,
  stopped before the last few) passed 64,725 cases, skipped 16, and
  failed 2. `exsrch` seed 101 is a palindromic query (`CG`) with
  `--strand both`. vsearch orders its plus and minus hits on one target
  by `std::sort`. That sort is not stable. The order is therefore
  libc++'s. `clust` seed 95 is an empty sequence. vsearch's
  `--cluster_fast` CIGAR for it is `=`, and clust's is `339I`. The
  `uc-fast` comparator now matches them.
