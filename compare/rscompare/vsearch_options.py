"""vsearch's options and a parser for its command lines.

The option names and whether each takes a value come from the option
table in vsearch's `src/cli.cc` (`../vsearch`, e1551256).
"""

from typing import Dict, List, Optional

# options that take no value
FLAGS = set("""
acceptall allow_fewer bzip2_decompress centroid_sizeout clusterout_id
clusterout_sort cons_truncate eeout fasta_score fastq_allowmergestagger
fastq_eeout fastq_nostagger fastq_qout_max fastq_solexa fulldp
gzip_decompress h hardmask help label_substr_match leftjust lengthout
n_mismatch no_progress notrunclabels output_no_hits quiet relabel_keep
relabel_md5 relabel_self relabel_sha1 rightjust samheader self selfid
sff_clip sintax_random sizein sizeorder sizeout top_hits_only uc_allhits
uchimeout5 usersort v version xee xlength xsize
""".split())

# options that take a value, the commands included
VALUED = set("""
abskew alignwidth allpairs_global alnout band biomout blast6out borderline
centroids chimeras chimeras_denovo chimeras_diff_pct chimeras_length_min
chimeras_parents_max chimeras_parts cluster_fast cluster_size
cluster_smallmem cluster_unoise clusters consout cut cut_pattern db dbmask
dbmatched dbnotmatched derep_fulllength derep_id derep_prefix
derep_smallmem dn ee_cutoffs eetabbedout fasta2fastq fasta_width fastaout
fastaout_discarded fastaout_discarded_rev fastaout_notmerged_fwd
fastaout_notmerged_rev fastaout_orphans fastaout_orphans_rev fastaout_rev
fastapairs fastq_ascii fastq_asciiout fastq_chars fastq_convert
fastq_eestats fastq_eestats2 fastq_filter fastq_join fastq_maxdiffpct
fastq_maxdiffs fastq_maxee fastq_maxee_rate fastq_maxlen fastq_maxmergelen
fastq_maxns fastq_mergepairs fastq_minlen fastq_minmergelen fastq_minovlen
fastq_minqual fastq_qmax fastq_qmaxout fastq_qmin fastq_qminout fastq_stats
fastq_stripleft fastq_stripright fastq_tail fastq_truncee
fastq_truncee_rate fastq_trunclen fastq_trunclen_keep fastq_truncqual
fastqout fastqout_discarded fastqout_discarded_rev fastqout_notmerged_fwd
fastqout_notmerged_rev fastqout_orphans fastqout_orphans_rev fastqout_rev
fastx_filter fastx_getseq fastx_getseqs fastx_getsubseq fastx_mask
fastx_revcomp fastx_subsample fastx_syncpairs fastx_uniques gapext gapopen
hspw id iddef idprefix idsuffix join_padgap join_padgapq label label_field
label_suffix label_word label_words labels lca_cutoff lcaout length_cutoffs
log makeudb_usearch maskfasta match matched max_unmasked_pct maxaccepts
maxdiffs maxgaps maxhits maxid maxqsize maxqt maxrejects maxseqlength
maxsize maxsizeratio maxsl maxsubs maxuniquesize mid min_unmasked_pct
mincols mindiffs mindiv minh minhsp minqt minseqlength minsize
minsizeratio minsl mintsize minuniquesize minwordmatches mismatch
mothur_shared_out msaout nonchimeras notmatched notmatchedfq orient
otutabout output pattern profile qmask qsegout query_cov randseed
read_separators relabel rereplicate reverse rowlen samout sample
sample_pct sample_size scramble scramble_kmer search_exact sff_convert
shuffle sintax sintax_cutoff slots sortbylength sortbysize strand
subseq_end subseq_start tabbedout target_cov threads topn tsegout uc
uchime2_denovo uchime3_denovo uchime_denovo uchime_ref uchimealns
uchimeout udb2fasta udbinfo udbstats unoise_alpha usearch_global
userfields userout weak_id wordlength xdrop_nw xn
""".split())

# the 48 commands (and --help, --version, -h, -v)
COMMANDS = set("""
allpairs_global chimeras_denovo cluster_fast cluster_size cluster_smallmem
cluster_unoise cut derep_fulllength derep_id derep_prefix derep_smallmem
fasta2fastq fastq_chars fastq_convert fastq_eestats fastq_eestats2
fastq_filter fastq_join fastq_mergepairs fastq_stats fastx_filter
fastx_getseq fastx_getseqs fastx_getsubseq fastx_mask fastx_revcomp
fastx_subsample fastx_syncpairs fastx_uniques makeudb_usearch maskfasta
orient rereplicate scramble search_exact sff_convert shuffle sintax sortbylength
sortbysize uchime2_denovo uchime3_denovo uchime_denovo uchime_ref
udb2fasta udbinfo udbstats usearch_global
help version h v
""".split())


class TranslationError(Exception):
    """A vsearch command line the translator cannot turn into rsearch."""


class Args:
    """A parsed vsearch command line: the command, its argument, and the
    other options in order. Options are taken by the translator as it uses
    them; whatever is left untaken is untranslatable."""

    def __init__(self, argv: List[str]):
        if argv and not argv[0].startswith('-'):
            argv = argv[1:]  # the program name
        self.argv = list(argv)
        self.opts: Dict[str, object] = {}
        self.order: List[str] = []
        self.command: Optional[str] = None
        self.input: Optional[str] = None
        self.used = set()
        i = 0
        while i < len(argv):
            tok = argv[i]
            if not tok.startswith('-') or tok in ('-', '--'):
                raise TranslationError(f'unexpected argument {tok!r}')
            name = tok.lstrip('-')
            value = None
            if '=' in name:
                name, value = name.split('=', 1)
            if name in FLAGS and name not in COMMANDS:
                if value is not None:
                    raise TranslationError(f'--{name} takes no value')
                value = True
            elif name in VALUED or name in COMMANDS:
                if value is None and name not in ('help', 'version', 'h', 'v'):
                    if i + 1 >= len(argv):
                        raise TranslationError(f'--{name} needs a value')
                    i += 1
                    value = argv[i]
                elif value is None:
                    value = True
            else:
                raise TranslationError(f'unknown vsearch option --{name}')
            if name in COMMANDS:
                if self.command is not None:
                    raise TranslationError(
                        f'two commands: --{self.command} and --{name}')
                self.command = name
                self.input = value if value is not True else None
            else:
                if name not in self.opts:
                    self.order.append(name)
                self.opts[name] = value
            i += 1
        if self.command is None:
            raise TranslationError('no vsearch command')

    def has(self, name: str) -> bool:
        return name in self.opts

    def take(self, name: str, default=None):
        """The option's value (marking it used), or the default."""
        if name in self.opts:
            self.used.add(name)
            return self.opts[name]
        return default

    def flag(self, name: str) -> bool:
        return bool(self.take(name, False))

    def take_num(self, name: str, default=None):
        v = self.take(name)
        if v is None:
            return default
        try:
            return float(v) if any(c in str(v) for c in '.eE') else int(v)
        except ValueError:
            raise TranslationError(f'--{name} {v}: not a number')

    def untaken(self) -> List[str]:
        return [n for n in self.order if n not in self.used]
