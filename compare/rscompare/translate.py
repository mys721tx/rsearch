"""Translates a vsearch command line into the rsearch pipeline that gives
the same results.

The rules come from the pipelines of `README.md`, the tool pages in
`docs/` and the option map in `docs/vsearch.md`, each written once here:

    >>> p = translate(['--fastx_filter', 'in.fq', '--fastq_maxee', '1',
    ...                '--fastqout', 'out.fq'], 'work')
    >>> print(p.render())
    sqfilt -i in.fq --max-ee 1 > work/fastqout

`translate(argv, workdir)` returns a `Pipeline`: the rsearch commands in
the order they run (`steps`, each a tool, its arguments and the files on
its standard input and output), and for each output of the vsearch
command the rsearch file that stands for it and the comparator that
compares them (`outputs`). Temporary files and rsearch's outputs go into
`workdir`. Steps run one after the other, each reading the files the
earlier ones wrote, which is a shell pipeline without the concurrency.

An option it cannot translate raises `TranslationError` naming it. A
behaviour that is vsearch's default rather than an option, and that the
pipeline does not reproduce (a length filter sqfilt cannot express, labels
cut at spaces without `--notrunclabels`, several threads), is listed in
`caveats` instead.
"""

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .vsearch_options import Args, TranslationError

__all__ = ['translate', 'Pipeline', 'Step', 'Output', 'TranslationError']


@dataclass
class Step:
    tool: str
    args: List[str]
    stdin: Optional[str] = None
    stdout: Optional[str] = None

    def argv(self, bindir: str) -> List[str]:
        return [os.path.join(bindir, self.tool)] + self.args

    def inputs(self) -> List[str]:
        """The files this step reads (its standard input and -i)."""
        out = [self.stdin] if self.stdin else []
        for k, a in enumerate(self.args[:-1]):
            if a in ('-i', '--reverse'):
                out.append(self.args[k + 1])
        return out


@dataclass
class Output:
    option: str        # the vsearch option, without dashes ("fastaout")
    vsearch: str       # its file (or prefix, for --clusters) in vsearch's argv
    rsearch: str       # the rsearch file (or prefix) that stands for it
    comparator: str = 'bytes'


@dataclass
class Pipeline:
    command: str
    argv: List[str]
    steps: List[Step] = field(default_factory=list)
    outputs: List[Output] = field(default_factory=list)
    caveats: List[str] = field(default_factory=list)

    def render(self) -> str:
        """The steps as shell pipelines, a line each."""
        lines = []
        uses: Dict[str, int] = {}
        for s in self.steps:
            for f in s.inputs():
                uses[f] = uses.get(f, 0) + 1
        pending = None
        for s in self.steps:
            words = [s.tool] + [_q(a) for a in s.args]
            if pending is not None and s.stdin == pending[1] \
                    and uses.get(s.stdin, 0) == 1 and _tmp(s.stdin):
                text = pending[0] + ' | ' + ' '.join(words)
            else:
                if pending is not None:
                    lines.append(_close(pending))
                text = ' '.join(words)
                if s.stdin:
                    text += ' < ' + s.stdin
            pending = (text, s.stdout)
        if pending is not None:
            lines.append(_close(pending))
        return '\n'.join(lines)

    def to_json(self) -> dict:
        return {'command': self.command,
                'steps': [{'tool': s.tool, 'args': s.args, 'stdin': s.stdin,
                           'stdout': s.stdout} for s in self.steps],
                'outputs': [o.__dict__ for o in self.outputs],
                'caveats': self.caveats}


def _q(a: str) -> str:
    if a == '' or any(c in a for c in ' \t\'"$*;|&<>()^\\'):
        return "'" + a.replace("'", "'\\''") + "'"
    return a


def _tmp(path: Optional[str]) -> bool:
    return bool(path) and os.path.basename(path).startswith('tmp')


def _close(pending) -> str:
    text, out = pending
    return text + (' > ' + out if out else '')


# the first bytes of a file: FASTA, FASTQ, UDB or SFF, or None
def sniff(path: Optional[str]) -> Optional[str]:
    if not path or path == '-':
        return None
    try:
        with open(path, 'rb') as f:
            head = f.read(4)
    except OSError:
        return None
    if head[:1] == b'>':
        return 'fasta'
    if head[:1] == b'@':
        return 'fastq'
    if head == b'FBDU':
        return 'udb'
    if head == b'.sff':
        return 'sff'
    return None


# acceptance criteria and scores: vsearch name -> rsearch option
CRITERIA = {
    'id': '--id', 'iddef': '--iddef', 'weak_id': '--weak-id',
    'maxdiffs': '--max-diffs', 'maxsubs': '--max-subs',
    'maxgaps': '--max-gaps', 'mincols': '--min-cols',
    'query_cov': '--query-cov', 'target_cov': '--target-cov',
    'maxid': '--max-id', 'mid': '--mid', 'idprefix': '--id-prefix',
    'idsuffix': '--id-suffix', 'maxqt': '--max-qt', 'minqt': '--min-qt',
    'minsl': '--min-sl', 'maxsl': '--max-sl', 'maxqsize': '--max-qsize',
    'mintsize': '--min-tsize', 'minsizeratio': '--min-size-ratio',
    'maxsizeratio': '--max-size-ratio',
    'match': '--match', 'mismatch': '--mismatch',
    'gapopen': '--gapopen', 'gapext': '--gapext',
}
CRITERIA_FLAGS = {'leftjust': '--left-just', 'rightjust': '--right-just',
                  'self': '--self', 'selfid': '--self-id',
                  'n_mismatch': '--n-mismatch'}
WORDS = {'strand': '--strand', 'maxaccepts': '--maxaccepts',
         'maxrejects': '--maxrejects', 'wordlength': '--word-len',
         'minwordmatches': '--min-word-matches'}
# options vsearch accepts and ignores
IGNORED = ['band', 'hspw', 'minhsp', 'xdrop_nw', 'fulldp', 'slots', 'pattern']
# hitfmt takes the options the search scored with (decision 25)
HITFMT_SCORING = {'iddef': '--iddef', 'match': '--match',
                  'mismatch': '--mismatch', 'gapopen': '--gapopen',
                  'gapext': '--gapext'}

LABELS = {'relabel': '--relabel', 'sample': '--sample',
          'label_suffix': '--label-suffix'}
LABEL_FLAGS = {'relabel_keep': '--relabel-keep',
               'relabel_md5': '--relabel-md5',
               'relabel_self': '--relabel-self',
               'relabel_sha1': '--relabel-sha1', 'sizeout': '--sizeout',
               'xsize': '--xsize', 'lengthout': '--lengthout',
               'xlength': '--xlength', 'xee': '--xee', 'eeout': '--ee',
               'fastq_eeout': '--ee'}

THREADED = {'usearch_global', 'search_exact', 'allpairs_global',
            'cluster_fast', 'cluster_size', 'cluster_smallmem',
            'cluster_unoise', 'uchime_ref', 'sintax', 'fastq_mergepairs',
            'fastx_mask', 'maskfasta', 'makeudb_usearch'}
# commands that keep whole headers without --notrunclabels, or read none
WHOLE_HEADERS = {'fastq_convert', 'fasta2fastq', 'fastx_filter',
                 'fastq_filter', 'fastx_revcomp', 'sintax', 'fastq_chars',
                 'fastq_stats', 'fastq_eestats', 'fastq_eestats2',
                 'udbinfo', 'udbstats', 'udb2fasta', 'sff_convert',
                 'fastq_join', 'fastq_mergepairs', 'fastx_syncpairs',
                 'cut', 'rereplicate', 'scramble', 'fastx_getseq',
                 'fastx_getseqs', 'fastx_getsubseq'}


class _Translator:
    def __init__(self, argv: List[str], workdir: str):
        self.a = Args(argv)
        self.cmd = self.a.command
        self.inp = self.a.input
        self.work = workdir
        self.p = Pipeline(self.cmd, list(argv))
        self.n = 0
        self.width = 80

    # building blocks

    def caveat(self, text: str):
        if text not in self.p.caveats:
            self.p.caveats.append(text)

    def tmp(self, name: str) -> str:
        self.n += 1
        return os.path.join(self.work, f'tmp{self.n}.{name}')

    def path(self, name: str) -> str:
        return os.path.join(self.work, name)

    def run(self, tool: str, args: List[str], stdin: Optional[str] = None,
            stdout: Optional[str] = None, name: str = 'out') -> str:
        """Adds a step; returns the file its standard output goes to."""
        if stdout is None:
            stdout = self.tmp(name)
        self.p.steps.append(Step(tool, list(args), stdin, stdout))
        return stdout

    def chain(self, src: str, cmds: List[Tuple[str, List[str]]],
              first_reads_file: bool = True) -> str:
        """Runs the commands one into the next, from src."""
        for tool, args in cmds:
            src = self.step(tool, src, args, name=tool)
        return src

    def step(self, tool: str, src: str, args: List[str],
             stdout: Optional[str] = None, name: str = 'out') -> str:
        """A step reading src: on its standard input when src is a
        temporary file, else as -i."""
        if _tmp(src):
            return self.run(tool, args, stdin=src, stdout=stdout, name=name)
        return self.run(tool, ['-i', src] + args, stdout=stdout, name=name)

    def output(self, option: str, rpath: str, comparator: str = 'bytes'):
        self.p.outputs.append(Output(option, self.a.take(option), rpath,
                                     comparator))

    # outputs of records

    def fasta(self, option: str, src: str, fmt: Optional[str]):
        """vsearch's FASTA output: conv --fasta at --fasta_width, unless the
        stream is FASTA already and unwrapped."""
        if not self.a.has(option):
            return
        if fmt == 'fasta' and self.width == 0:
            self.output(option, src)
            return
        args = ['--fasta'] + ([] if self.width == 80 else
                              ['--fasta-width', str(self.width)])
        if _tmp(src):
            out = self.run('conv', args, stdin=src, name=option)
        else:
            out = self.run('conv', ['-i', src] + args, name=option)
        self.output(option, out)

    def fastq(self, option: str, src: str):
        if self.a.has(option):
            self.output(option, src)

    def records(self, src: str, fmt: Optional[str], fa: Optional[str],
                fq: Optional[str] = None):
        if fa:
            self.fasta(fa, src, fmt)
        if fq:
            self.fastq(fq, src)

    def unwrapped_only(self, option: str):
        """Outputs rsearch writes unwrapped, which conv cannot wrap."""
        if self.a.has(option) and self.width != 0:
            raise TranslationError(
                f'--{option} with --fasta_width {self.width}: rsearch writes '
                'it unwrapped (give --fasta_width 0)')

    # options shared by many commands

    def qual_in(self) -> List[str]:
        out = []
        for v, r in (('fastq_ascii', '--ascii'), ('fastq_qmin', '--qmin'),
                     ('fastq_qmax', '--qmax')):
            if self.a.has(v):
                out += [r, str(self.a.take(v))]
        return out

    def qual_out(self) -> List[str]:
        out = []
        for v, r in (('fastq_asciiout', '--ascii-out'),
                     ('fastq_qminout', '--qmin-out'),
                     ('fastq_qmaxout', '--qmax-out')):
            if self.a.has(v):
                out += [r, str(self.a.take(v))]
        return out

    def label_args(self, offset: Optional[str] = None,
                   skip: Tuple[str, ...] = ()) -> List[str]:
        out = []
        for v, r in LABELS.items():
            if v not in skip and self.a.has(v):
                out += [r, str(self.a.take(v))]
        for v, r in LABEL_FLAGS.items():
            if v not in skip and self.a.flag(v) and r not in out:
                out.append(r)
        if '--ee' in out and offset and offset != '33':
            out += ['--ascii', offset]
        return out

    def labels(self, src: str, offset: Optional[str] = None,
               skip: Tuple[str, ...] = ()) -> str:
        args = self.label_args(offset, skip)
        return self.run('relab', args, stdin=src, name='relab') if args else src

    def no_labels(self, why: str, skip: Tuple[str, ...] = ()):
        for v in list(LABELS) + list(LABEL_FLAGS):
            if v not in skip and self.a.has(v) and v not in self.a.used:
                raise TranslationError(f'--{v} with --{self.cmd}: {why}')

    def lengths(self, min_default: int) -> Tuple[int, int]:
        return (self.a.take_num('minseqlength', min_default),
                self.a.take_num('maxseqlength', 50000))

    def length_filter(self, min_default: int, extra: List[str] = (),
                      what: str = 'input') -> Optional[List[str]]:
        """sqfilt's arguments for vsearch's length filters (and any other
        filter of the same step), or None for no sqfilt."""
        mn, mx = self.lengths(min_default)
        extra = list(extra)
        if mn <= 0:
            if extra:
                raise TranslationError(
                    f'--minseqlength {mn} with {" ".join(extra)}: sqfilt '
                    'cannot keep empty sequences')
            self.caveat(f'--maxseqlength {mx} not applied to the {what}: '
                        f'with --minseqlength {mn} the empty sequences are '
                        'kept, which sqfilt cannot do')
            return None
        args = [] if mn == 1 else ['--min-len', str(mn)]
        return args + ['--max-len', str(mx)] + extra

    def filtered(self, src: str, min_default: int, extra: List[str] = (),
                 what: str = 'input') -> str:
        args = self.length_filter(min_default, extra, what)
        if args is None:
            return src
        return self.run('sqfilt', ['-i', src] + args, name='sqfilt')

    def masking(self, option: str, style: str) -> List[Tuple[str, List[str]]]:
        """lcmask runs for --qmask or --dbmask (and --hardmask), decision 8
        and the checks under the pipeline tables."""
        mode = self.a.take(option, 'dust')
        hard = self.a.flag('hardmask')
        if mode not in ('dust', 'soft', 'none'):
            raise TranslationError(f'--{option} {mode}: not dust, soft or none')
        if style == 'udb':  # --makeudb_usearch upper-cases first
            if mode == 'dust':
                return [('lcmask', ['--clear']), ('lcmask', ['--hard'])] \
                    if hard else [('lcmask', [])]
            return [('lcmask', ['--clear'])]
        if mode == 'dust':
            return [('lcmask', ['--hard'] if hard else [])]
        if mode == 'soft':
            return [('lcmask', ['--hard', '--no-dust'])] if hard else []
        return [('lcmask', ['--clear'])]

    def masked(self, src: str, cmds) -> str:
        return self.chain(src, cmds) if cmds else src

    def search_args(self, words: bool = True, hitchoice: bool = True,
                    skip: Tuple[str, ...] = ()) -> List[str]:
        out = []
        table = dict(CRITERIA)
        if words:
            table.update(WORDS)
        for v, r in table.items():
            if v not in skip and self.a.has(v):
                out += [r, str(self.a.take(v))]
        for v, r in CRITERIA_FLAGS.items():
            if v not in skip and self.a.flag(v):
                out.append(r)
        if hitchoice:
            if self.a.has('maxhits'):
                out += ['--maxhits', str(self.a.take('maxhits'))]
            if self.a.flag('top_hits_only'):
                out.append('--top-hits-only')
        for v in IGNORED:
            if self.a.has(v):
                self.a.take(v)
                self.caveat(f'--{v} ignored, as vsearch ignores it')
        return out

    def hitfmt_scoring(self) -> List[str]:
        out = []
        for v, r in HITFMT_SCORING.items():
            if self.a.has(v):
                out += [r, str(self.a.opts[v])]
        if self.a.has('n_mismatch'):
            out.append('--n-mismatch')
        return out

    def threads(self):
        t = self.a.take('threads')
        if self.cmd in THREADED and str(t) != '1':
            self.caveat('vsearch runs on several threads without --threads 1; '
                        'rsearch gives its single-threaded results')

    # the translation

    def translate(self) -> Pipeline:
        a = self.a
        for v in ('quiet', 'no_progress'):
            a.take(v)
        if not a.flag('notrunclabels') and self.cmd not in WHOLE_HEADERS:
            self.caveat('labels are whole (decision 9): vsearch cuts them at '
                        'the first space or tab without --notrunclabels')
        for v in ('gzip_decompress', 'bzip2_decompress'):
            if a.has(v):
                raise TranslationError(f'--{v}: rsearch reads no compressed '
                                       'input')
        self.width = a.take_num('fasta_width', 80)
        self.threads()
        handler = getattr(self, 't_' + self.cmd, None)
        if handler is None:
            raise TranslationError(f'--{self.cmd}: no rsearch translation')
        if self.inp is None:
            raise TranslationError(f'--{self.cmd} needs an input file')
        handler()
        if a.has('log') and 'log' not in a.used:
            a.take('log')  # a log of the run, not an output compared
        left = a.untaken()
        if left:
            raise TranslationError(
                ', '.join('--' + n for n in left) +
                f': no rsearch translation for --{self.cmd}')
        return self.p

    def log_output(self, tool: str, args: List[str], comparator: str,
                   stdin: Optional[str] = None):
        """A report that vsearch writes to --log: the tool's -o."""
        if self.a.has('log'):
            out = self.path('log')
            self.run(tool, args, stdin=stdin, stdout=out)
            self.output('log', out, comparator)
        else:
            self.run(tool, args, stdin=stdin, name=tool)

    # format conversion

    def t_fastq_convert(self):
        args = self.qual_in() + self.qual_out()
        if self.a.flag('fastq_solexa'):
            args.append('--solexa')
        offset = str(self.a.opts.get('fastq_asciiout',
                                     self.a.opts.get('fastq_ascii', 33)))
        s = self.run('conv', ['-i', self.inp] + args, name='conv')
        self.fastq('fastqout', self.labels(s, offset, skip=('sizein',)))
        self.a.take('sizein')

    def t_fasta2fastq(self):
        args = self.qual_out()
        offset = str(self.a.opts.get('fastq_asciiout', 33))
        s = self.run('conv', ['-i', self.inp] + args, name='conv')
        self.fastq('fastqout', self.labels(s, offset))
        self.a.take('sizein')

    def t_sff_convert(self):
        args = (['--clip'] if self.a.flag('sff_clip') else []) + self.qual_out()
        s = self.run('conv', ['-i', self.inp] + args, name='conv')
        self.fastq('fastqout', self.labels(s))

    def t_udb2fasta(self):
        if self.label_args():
            raise TranslationError('label options with --udb2fasta')
        self.fasta('output', self.inp, 'udb')

    # quality reports

    def t_fastq_chars(self):
        args = ['-i', self.inp]
        if self.a.has('fastq_tail'):
            args += ['--tail', str(self.a.take('fastq_tail'))]
        self.log_output('fqchar', args, 'log-report')

    def t_fastq_stats(self):
        self.log_output('fqstat', ['-i', self.inp] + self.qual_in(),
                        'log-report')

    def t_fastq_eestats(self):
        s = self.run('eestat', ['-i', self.inp] + self.qual_in(),
                     name='eestat')
        self.fastq('output', s)

    def t_fastq_eestats2(self):
        args = ['-i', self.inp]
        if self.a.has('length_cutoffs'):
            parts = str(self.a.take('length_cutoffs')).split(',')
            if len(parts) != 3:
                raise TranslationError('--length_cutoffs needs three values')
            args += ['--min-len', parts[0]]
            if parts[1] != '*':
                args += ['--max-len', parts[1]]
            args += ['--len-step', parts[2]]
        if self.a.has('ee_cutoffs'):
            args += ['--max-ee', str(self.a.take('ee_cutoffs'))]
        s = self.run('eegrid', args + self.qual_in(), name='eegrid')
        self.fastq('output', s)

    # trimming and filtering

    def t_fastx_filter(self):
        a = self.a
        trim = []
        for v, r in (('fastq_stripleft', '--strip-left'),
                     ('fastq_stripright', '--strip-right'),
                     ('fastq_trunclen_keep', '--len'),
                     ('fastq_truncqual', '--qual'),
                     ('fastq_truncee', '--max-ee'),
                     ('fastq_truncee_rate', '--max-ee-rate')):
            if a.has(v):
                trim += [r, str(a.take(v))]
        minlen = a.take_num('fastq_minlen')
        if a.has('fastq_trunclen'):
            if a.has('fastq_trunclen_keep'):
                raise TranslationError('--fastq_trunclen with '
                                       '--fastq_trunclen_keep')
            n = a.take_num('fastq_trunclen')
            trim += ['--len', str(n)]
            minlen = n if minlen is None else max(n, minlen)
        filt = []
        if minlen is not None:
            if minlen < 1:
                raise TranslationError(f'--fastq_minlen {minlen}: sqfilt '
                                       'keeps no empty read')
            filt += ['--min-len', str(minlen)]
        for v, r in (('fastq_maxlen', '--max-len'), ('fastq_maxns', '--max-ns'),
                     ('fastq_maxee', '--max-ee'),
                     ('fastq_maxee_rate', '--max-ee-rate'),
                     ('fastq_minqual', '--min-qual'), ('minsize', '--min-size'),
                     ('maxsize', '--max-size')):
            if a.has(v):
                filt += [r, str(a.take(v))]
        offset = str(a.opts.get('fastq_ascii', 33))
        q = self.qual_in()
        rev = a.take('reverse')
        fmt = sniff(self.inp)
        src, rsrc = self.inp, rev
        if trim:
            src = self.run('sqtrim', ['-i', src] + trim + q, name='sqtrim')
            if rev:
                rsrc = self.run('sqtrim', ['-i', rsrc] + trim + q,
                                name='sqtrim')
        kept = self.tmp('kept')
        args = []
        disc = rkept = rdisc = None
        if rev:
            rkept = self.tmp('rev_kept')
            args += ['--reverse', rsrc, '--rev-out', rkept]
        wants = lambda *o: any(a.has(x) for x in o)
        if wants('fastaout_discarded', 'fastqout_discarded'):
            disc = self.tmp('discarded')
            args += ['--discarded', disc]
        if rev and wants('fastaout_discarded_rev', 'fastqout_discarded_rev'):
            rdisc = self.tmp('rev_discarded')
            args += ['--rev-discarded', rdisc]
        self.step('sqfilt', src, args + filt + q, stdout=kept)
        a.take('sizein')
        if disc or rev:
            self.no_labels('relabelling with several outputs is not '
                           'translated')
        kept = self.labels(kept, offset)
        self.records(kept, fmt, 'fastaout', 'fastqout')
        if disc:
            self.records(disc, fmt, 'fastaout_discarded', 'fastqout_discarded')
        if rev:
            rfmt = sniff(rev)
            self.records(rkept, rfmt, 'fastaout_rev', 'fastqout_rev')
            if rdisc:
                self.records(rdisc, rfmt, 'fastaout_discarded_rev',
                             'fastqout_discarded_rev')

    t_fastq_filter = t_fastx_filter

    # pairs

    def t_fastq_mergepairs(self):
        a = self.a
        rev = a.take('reverse')
        if rev is None:
            raise TranslationError('--fastq_mergepairs needs --reverse')
        q = self.qual_in()
        offset = str(a.opts.get('fastq_ascii', 33))
        fwd = self.inp
        pre = []
        if a.has('fastq_truncqual'):
            t = ['--qual', str(a.take('fastq_truncqual'))]
            fwd = self.run('sqtrim', ['-i', fwd] + t + q, name='sqtrim')
            rev = self.run('sqtrim', ['-i', rev] + t + q, name='sqtrim')
            pre.append('--fastq_truncqual')
        filt = []
        for v, r in (('fastq_minlen', '--min-len'), ('fastq_maxlen', '--max-len'),
                     ('fastq_maxns', '--max-ns')):
            if a.has(v):
                filt += [r, str(a.take(v))]
                pre.append('--' + v)
        if filt:
            rkept = self.tmp('rev')
            fwd = self.run('sqfilt', ['-i', fwd, '--reverse', rev,
                                      '--rev-out', rkept] + filt + q,
                           name='sqfilt')
            rev = rkept
        margs = ['-i', fwd, '--reverse', rev]
        for v, r in (('fastq_minovlen', '--min-overlap'),
                     ('fastq_maxdiffs', '--max-diffs'),
                     ('fastq_maxdiffpct', '--max-diff-pct'),
                     ('fastq_qminout', '--qmin-out'),
                     ('fastq_qmaxout', '--qmax-out')):
            if a.has(v):
                margs += [r, str(a.take(v))]
        if a.flag('fastq_allowmergestagger'):
            margs.append('--stagger')
        a.take('fastq_nostagger')
        post = []
        for v, r in (('fastq_maxee', '--max-ee'),
                     ('fastq_minmergelen', '--min-len'),
                     ('fastq_maxmergelen', '--max-len')):
            if a.has(v):
                post += [r, str(a.take(v))]
        filters = pre + [x for x in post if x.startswith('--')]
        side = {}
        for opt, flag in (('notmerged_fwd', '--not-merged'),
                          ('notmerged_rev', '--rev-not-merged')):
            if a.has('fastaout_' + opt) or a.has('fastqout_' + opt):
                if filters:
                    raise TranslationError(
                        f'--fast?out_{opt} with {", ".join(filters)}: the '
                        'reads sqfilt discards are not written as not merged')
                side[opt] = self.tmp(opt)
                margs += [flag, side[opt]]
        if a.has('eetabbedout'):
            if post:
                raise TranslationError('--eetabbedout with filters of merged '
                                       'reads')
            margs += ['--ee-table', self.path('eetabbedout')]
        report = None
        if a.has('log'):
            if filters:
                a.take('log')
                self.caveat('the --log report is not compared: filters run '
                            'outside mergpa')
            else:
                report = self.path('log')
                margs += ['--report', report]
        merged = self.run('mergpa', margs + q, name='merged')
        if a.has('eetabbedout'):
            self.output('eetabbedout', self.path('eetabbedout'))
        if report:
            self.output('log', report, 'log-report')
        if post:
            merged = self.run('sqfilt', post + q, stdin=merged, name='sqfilt')
        a.take('sizein')
        merged = self.labels(merged, offset)
        self.records(merged, 'fastq', 'fastaout', 'fastqout')
        for opt, f in side.items():
            self.records(f, 'fastq', 'fastaout_' + opt, 'fastqout_' + opt)

    def t_fastq_join(self):
        a = self.a
        rev = a.take('reverse')
        if rev is None:
            raise TranslationError('--fastq_join needs --reverse')
        args = ['-i', self.inp, '--reverse', rev]
        if a.has('join_padgap'):
            args += ['--pad', str(a.take('join_padgap'))]
        if a.has('join_padgapq'):
            args += ['--pad-qual', str(a.take('join_padgapq'))]
        if a.has('fastq_ascii'):
            args += ['--ascii', str(a.take('fastq_ascii'))]
        s = self.run('joinpa', args, name='joinpa')
        a.take('sizein')
        s = self.labels(s, str(a.opts.get('fastq_ascii', 33)))
        self.records(s, 'fastq', 'fastaout', 'fastqout')

    def t_fastx_syncpairs(self):
        a = self.a
        rev = a.take('reverse')
        if rev is None:
            raise TranslationError('--fastx_syncpairs needs --reverse')
        outs = {k: self.tmp(k) for k in ('fwd', 'rev', 'orphans', 'rev_orphans')}
        args = ['-i', self.inp, '--reverse', rev, '-o', outs['fwd'],
                '--rev-out', outs['rev'], '--orphans', outs['orphans'],
                '--rev-orphans', outs['rev_orphans']]
        if a.has('read_separators'):
            args += ['--separators', str(a.take('read_separators'))]
        self.run('syncpa', args, name='syncpa')
        ffmt, rfmt = sniff(self.inp), sniff(rev)
        self.records(outs['fwd'], ffmt, 'fastaout', 'fastqout')
        self.records(outs['rev'], rfmt, 'fastaout_rev', 'fastqout_rev')
        self.records(outs['orphans'], ffmt, 'fastaout_orphans',
                     'fastqout_orphans')
        self.records(outs['rev_orphans'], rfmt, 'fastaout_orphans_rev',
                     'fastqout_orphans_rev')

    # dereplication

    def derep(self, mode: List[str], min_default: int, fa: str,
              fq: Optional[str] = None):
        a = self.a
        fmt = sniff(self.inp)
        src = self.filtered(self.inp, min_default)
        args = list(mode)
        if a.has('strand'):
            s = a.take('strand')
            if s == 'both':
                args.append('--strand')
            elif s != 'plus':
                raise TranslationError(f'--strand {s}')
        sizein, sizeout = a.flag('sizein'), a.flag('sizeout')
        if not sizein:
            args.append('--no-sizein')
            if src != self.inp:
                self.caveat('sqfilt checks ;size= annotations, which '
                            'dereplication reads only with --sizein')
        if not sizeout:
            args.append('--no-sizeout')
        for v, r in (('uc', '--uc'), ('tabbedout', '--tabbedout')):
            if a.has(v):
                args += [r, self.path(v)]
        args += self.qual_in() + self.qual_out()
        if a.flag('fastq_qout_max'):
            args.append('--best-qual')
        s = self.step('derep', src, args, name='derep')
        for v in ('uc', 'tabbedout'):
            if a.has(v):
                self.output(v, self.path(v))
        sizes = []
        for v, r in (('minuniquesize', '--min-size'),
                     ('maxuniquesize', '--max-size')):
            if a.has(v):
                sizes += [r, str(a.take(v))]
        if sizes:
            if not sizeout:
                raise TranslationError('--minuniquesize or --maxuniquesize '
                                       'without --sizeout')
            s = self.run('sqfilt', sizes, stdin=s, name='sqfilt')
        if a.has('topn'):
            s = self.run('topn', [str(a.take('topn'))], stdin=s, name='topn')
        offset = str(a.opts.get('fastq_asciiout', a.opts.get('fastq_ascii', 33)))
        s = self.labels(s, offset, skip=('sizeout',))
        self.records(s, fmt, fa, fq)

    def t_fastx_uniques(self):
        self.derep([], 1, 'fastaout', 'fastqout')

    def t_derep_fulllength(self):
        self.derep([], 32, 'output')

    def t_derep_id(self):
        self.derep(['--same-label'], 32, 'output')

    def t_derep_prefix(self):
        self.derep(['--prefix'], 32, 'output')

    def t_derep_smallmem(self):
        self.derep(['--unsorted'], 1, 'fastaout')

    def t_rereplicate(self):
        a = self.a
        args = ['-i', self.inp] + ([] if a.flag('sizeout') else ['--no-sizeout'])
        a.take('sizein')
        s = self.run('rerep', args, name='rerep')
        s = self.labels(s, skip=('sizeout',))
        self.fasta('output', s, sniff(self.inp))

    # sequence manipulation

    def t_fastx_revcomp(self):
        s = self.run('revcomp', ['-i', self.inp], name='revcomp')
        self.a.take('sizein')
        s = self.labels(s, str(self.a.opts.get('fastq_ascii', 33)))
        self.records(s, sniff(self.inp), 'fastaout', 'fastqout')

    def mask_command(self, fa: str, fq: Optional[str]):
        a = self.a
        mode = a.take('qmask', 'dust')
        hard = a.flag('hardmask')
        if mode == 'dust':
            cmds = [('lcmask', ['--hard'] if hard else [])]
        elif mode == 'soft':
            cmds = [('lcmask', ['--hard', '--no-dust'])] if hard else []
        elif mode == 'none':
            cmds = []
        else:
            raise TranslationError(f'--qmask {mode}')
        pct = []
        for v, r in (('min_unmasked_pct', '--min-%s-pct'),
                     ('max_unmasked_pct', '--max-%s-pct')):
            if a.has(v):
                if mode == 'none':
                    raise TranslationError(f'--{v} with --qmask none')
                pct += [r % ('non-n' if hard else 'upper'), str(a.take(v))]
        s = self.masked(self.inp, cmds)
        args = self.length_filter(1, pct)
        if args is not None:
            s = self.run('sqfilt', args, stdin=s, name='sqfilt') if cmds else \
                self.run('sqfilt', ['-i', s] + args, name='sqfilt')
        a.take('sizein')
        s = self.labels(s, str(a.opts.get('fastq_ascii', 33)))
        self.records(s, sniff(self.inp), fa, fq)

    def t_fastx_mask(self):
        self.mask_command('fastaout', 'fastqout')

    def t_maskfasta(self):
        self.mask_command('output', None)

    def t_cut(self):
        a = self.a
        if not a.has('cut_pattern'):
            raise TranslationError('--cut needs --cut_pattern')
        self.no_labels('relabelling the four outputs is not translated',
                       skip=('sizein',))
        a.take('sizein')
        outs = {k: self.tmp(k) for k in ('cut', 'rev', 'uncut', 'rev_uncut')}
        self.run('redig', ['-i', self.inp, '--pattern',
                           str(a.take('cut_pattern')), '--rev-out', outs['rev'],
                           '--uncut', outs['uncut'], '--rev-uncut',
                           outs['rev_uncut']], stdout=outs['cut'])
        fmt = sniff(self.inp)
        self.fasta('fastaout', outs['cut'], fmt)
        self.fasta('fastaout_rev', outs['rev'], fmt)
        self.fasta('fastaout_discarded', outs['uncut'], fmt)
        self.fasta('fastaout_discarded_rev', outs['rev_uncut'], fmt)

    def seed(self) -> List[str]:
        s = self.a.take('randseed')
        if s is None or str(s) == '0':
            self.caveat('no --randseed: vsearch and rsearch draw different '
                        'random numbers')
            return []
        return ['--seed=' + str(s)]

    def t_scramble(self):
        args = ['-i', self.inp] + self.seed()
        if self.a.has('scramble_kmer'):
            args += ['--kmer', str(self.a.take('scramble_kmer'))]
        s = self.run('scrmbl', args, name='scrmbl')
        self.a.take('sizein')
        s = self.labels(s, str(self.a.opts.get('fastq_ascii', 33)))
        self.records(s, sniff(self.inp), 'fastaout', 'fastqout')

    # retrieval

    def select_args(self) -> List[str]:
        a = self.a
        args = []
        for v, r in (('label', '--label'), ('labels', '--labels'),
                     ('label_word', '--word'), ('label_words', '--words'),
                     ('label_field', '--field')):
            if a.has(v):
                args += [r, str(a.take(v))]
        if a.flag('label_substr_match'):
            args.append('--substr')
        return args

    def getseq(self, region: bool):
        a = self.a
        sel = self.select_args()
        if region and '--label' not in sel:
            raise TranslationError('--fastx_getsubseq needs --label')
        kept = self.tmp('kept')
        args = ['-i', self.inp] + sel
        disc = None
        if a.has('notmatched') or a.has('notmatchedfq'):
            disc = self.tmp('notmatched')
            args += ['--discarded', disc]
            self.no_labels('relabelling with --notmatched is not translated',
                           skip=('sizein',))
        self.run('getseq', args, stdout=kept)
        s = kept
        if region:
            sub = []
            for v, r in (('subseq_start', '--start'), ('subseq_end', '--end')):
                if a.has(v):
                    sub += [r, str(a.take(v))]
            s = self.run('subseq', sub, stdin=s, name='subseq')
        a.take('sizein')
        fmt = sniff(self.inp)
        s = self.labels(s, str(a.opts.get('fastq_ascii', 33)))
        self.records(s, fmt, 'fastaout', 'fastqout')
        if disc:
            self.fasta('notmatched', disc, fmt)
            self.fastq('notmatchedfq', disc)

    def t_fastx_getseq(self):
        self.getseq(False)

    t_fastx_getseqs = t_fastx_getseq

    def t_fastx_getsubseq(self):
        self.getseq(True)

    # sorting and sampling

    def sort(self, key: str):
        a = self.a
        extra = []
        if key == '--size':
            for v, r in (('minsize', '--min-size'), ('maxsize', '--max-size')):
                if a.has(v):
                    extra += [r, str(a.take(v))]
        s = self.filtered(self.inp, 1, extra)
        s = self.chain(s, [('sortby', [key])], first_reads_file=True)
        if a.has('topn'):
            s = self.run('topn', [str(a.take('topn'))], stdin=s, name='topn')
        a.take('sizein')
        s = self.labels(s)
        self.fasta('output', s, sniff(self.inp))

    def t_sortbylength(self):
        self.sort('--length')

    def t_sortbysize(self):
        self.sort('--size')

    def t_shuffle(self):
        a = self.a
        s = self.filtered(self.inp, 1)
        s = self.step('shufl', s, self.seed(), name='shufl')
        if a.has('topn'):
            s = self.run('topn', [str(a.take('topn'))], stdin=s, name='topn')
        a.take('sizein')
        s = self.labels(s)
        self.fasta('output', s, sniff(self.inp))

    def t_fastx_subsample(self):
        a = self.a
        s = self.filtered(self.inp, 1)
        args = []
        if a.has('sample_size'):
            args += ['--count', str(a.take('sample_size'))]
        if a.has('sample_pct'):
            args += ['--pct', str(a.take('sample_pct'))]
        sizein, sizeout = a.flag('sizein'), a.flag('sizeout')
        if sizein != sizeout:
            raise TranslationError('--sizein without --sizeout, or the other '
                                   'way round: subsmp --by-size does both')
        if sizein:
            args.append('--by-size')
        if a.flag('allow_fewer'):
            args.append('--allow-fewer')
        args += self.seed()
        disc = None
        if a.has('fastaout_discarded') or a.has('fastqout_discarded'):
            disc = self.tmp('discarded')
            args += ['--discarded', disc]
        self.no_labels('relabelling subsamples is not translated')
        kept = self.step('subsmp', s, args, name='subsmp')
        fmt = sniff(self.inp)
        self.records(kept, fmt, 'fastaout', 'fastqout')
        if disc:
            self.records(disc, fmt, 'fastaout_discarded', 'fastqout_discarded')

    # searches: hits and what is made from them

    def hit_outputs(self, hits, query: str, db: str, aligned: bool = True):
        """hitfmt, getseq, otutab and lcatax on the hits. `hits(kind)` gives
        the hit file: 'all' every hit the search keeps, 'top' its first hit
        a query (vsearch's --uc without --uc_allhits), 'every' without
        --maxhits or --top_hits_only (for --dbmatched, decision 28)."""
        a = self.a
        scoring = self.hitfmt_scoring()
        nohits = ['--output-no-hits'] if a.flag('output_no_hits') else []
        if a.has('uc'):
            self.output('uc', hits('all' if a.flag('uc_allhits') else 'top'))
        a.take('uc_allhits')
        fmts = [('alnout', ['--alignments'], 'alnout'),
                ('blast6out', ['--blast6'], 'bytes'),
                ('userout', None, 'bytes'),
                ('samout', ['--sam'], 'sam'),
                ('fastapairs', ['--fasta-pairs'], 'bytes'),
                ('qsegout', ['--qsegs'], 'bytes'),
                ('tsegout', ['--tsegs'], 'bytes')]
        for opt, args, comp in fmts:
            if not a.has(opt):
                continue
            if not aligned:
                raise TranslationError(f'--{opt}: no alignments to write')
            if opt == 'userout':
                if not a.has('userfields'):
                    raise TranslationError('--userout without --userfields')
                args = ['--userout', str(a.take('userfields'))]
            if opt == 'alnout' and a.has('rowlen'):
                args = args + ['--rowlen', str(a.take('rowlen'))]
            if opt == 'samout' and a.flag('samheader'):
                args = args + ['--sam-header']
            if opt == 'fastapairs':
                self.unwrapped_only(opt)
            out = self.run('hitfmt', [hits('all'), '-i', query, '--db', db]
                           + args + scoring + nohits, name=opt)
            if opt in ('qsegout', 'tsegout'):
                self.fasta(opt, out, 'fasta')
            else:
                self.output(opt, out, comp)
        for opt, kind, extra, src in (
                ('matched', 'all', [], query),
                ('notmatched', 'all', ['--invert'], query),
                ('dbmatched', 'every', ['--targets'], db),
                ('dbnotmatched', 'every', ['--targets', '--invert'], db)):
            if a.has(opt):
                if opt.startswith('db') and a.has('weak_id'):
                    self.caveat('weak hits count for getseq --targets, which '
                                'vsearch does not count for --dbmatched '
                                '(decision 28)')
                s = self.run('getseq', ['--hits', hits(kind), '-i', src] + extra,
                             name=opt)
                self.fasta(opt, s, sniff(src) if src != db else 'fasta')
        for opt, args, comp in (('otutabout', [], 'bytes'),
                                ('mothur_shared_out', ['--mothur'], 'bytes'),
                                ('biomout', ['--biom'], 'biom')):
            if a.has(opt):
                self.output(opt, self.run('otutab', [hits('all'), '--db', db]
                                          + args, name=opt), comp)
        if a.has('lcaout'):
            args = [hits('all'), '--db', db]
            if a.has('lca_cutoff'):
                args += ['--cutoff', str(a.take('lca_cutoff'))]
            self.output('lcaout', self.run('lcatax', args, name='lcaout'))

    def hit_files(self, tool: str, queries: str, base: List[str],
                  choice: List[str]):
        """The searches, queries streamed on stdin when made by a step."""
        cache = {}

        def hits(kind: str) -> str:
            key = {'all': tuple(choice), 'top': ('--maxhits', '1'),
                   'every': ()}[kind]
            if key not in cache:
                cache[key] = self.step(tool, queries, base + list(key),
                                       name='uc')
            return cache[key]
        return hits

    def no_search_labels(self):
        self.a.take('sizein')
        self.no_labels(f'--{self.cmd} relabels nothing rsearch writes')

    def database(self, min_default: int, mask_style: Optional[str],
                 what: str = 'database') -> Tuple[str, str]:
        """The database as searched and as written (hitfmt, getseq): sqfilt
        for the length filters, lcmask for --dbmask. A UDB file is taken as
        it is."""
        db = self.a.take('db')
        if db is None:
            raise TranslationError(f'--{self.cmd} needs --db')
        if sniff(db) == 'udb':
            self.lengths(min_default)
            if mask_style:
                self.a.take('dbmask')
                self.a.take('hardmask')
            return db, db
        mn, mx = self.lengths(min_default)
        filt = None if mn <= 0 else \
            ([] if mn == 1 else ['--min-len', str(mn)]) + ['--max-len', str(mx)]
        if mn <= 0:
            self.caveat(f'--maxseqlength {mx} not applied to the {what}: with '
                        f'--minseqlength {mn} the empty sequences are kept, '
                        'which sqfilt cannot do')
        mode = self.a.opts.get('dbmask', 'dust') if mask_style else None
        cmds = self.masking('dbmask', mask_style) if mask_style else []

        def build(cmds):
            s = self.masked(db, cmds)
            if filt is not None:
                s = self.run('sqfilt', (['-i', s] if s == db else []) + filt,
                             stdin=None if s == db else s, name='db')
            return s
        searched = build(cmds)
        shown = searched if mode != 'none' else build([])
        return searched, shown

    def queries(self, style: str = 'standard') -> Tuple[str, str]:
        mode = self.a.opts.get('qmask', 'dust')
        cmds = self.masking('qmask', style)
        searched = self.masked(self.inp, cmds)
        return searched, (self.inp if mode == 'none' else searched)

    def t_usearch_global(self):
        a = self.a
        if not a.has('id'):
            raise TranslationError('--usearch_global needs --id')
        q, qshown = self.queries()
        db, dbshown = self.database(32, 'standard')
        choice = []
        if a.has('maxhits'):
            choice += ['--maxhits', str(a.take('maxhits'))]
        if a.flag('top_hits_only'):
            choice.append('--top-hits-only')
        args = self.search_args(hitchoice=False)
        hits = self.hit_files('dbsrch', q, ['--db', db] + args, choice)
        self.no_search_labels()
        self.hit_outputs(hits, qshown, dbshown)

    def t_search_exact(self):
        a = self.a
        if a.flag('hardmask'):
            raise TranslationError('--hardmask with --search_exact')
        mode = a.take('qmask', 'dust')
        a.take('dbmask')
        db, _ = self.database(1, None)
        dbshown = db
        qshown = self.inp
        shows = any(a.has(o) for o in ('alnout', 'blast6out', 'userout',
                                       'samout', 'fastapairs', 'qsegout',
                                       'tsegout', 'matched', 'notmatched',
                                       'dbmatched', 'dbnotmatched'))
        if shows and (mode == 'dust' or a.opts.get('dbmask', 'dust') == 'dust'):
            self.caveat('the sequences --search_exact writes as vsearch '
                        'masked them: unverified')
            if mode == 'dust':
                qshown = self.chain(self.inp, [('lcmask', [])])
            if a.opts.get('dbmask', 'dust') == 'dust':
                dbshown = self.chain(db, [('lcmask', [])])
        base = ['--db', db]
        if a.has('strand'):
            base += ['--strand', str(a.take('strand'))]
        if a.flag('self'):
            base.append('--self')
        for v, r in (('maxqsize', '--max-qsize'), ('mintsize', '--min-tsize'),
                     ('minsizeratio', '--min-size-ratio'),
                     ('maxsizeratio', '--max-size-ratio')):
            if a.has(v):
                base += [r, str(a.take(v))]
        choice = ['--maxhits', str(a.take('maxhits'))] \
            if a.has('maxhits') and str(a.opts['maxhits']) != '0' else []
        a.take('maxhits')
        a.take('top_hits_only')  # every exact hit is a top hit
        hits = self.hit_files('exsrch', self.inp, base, choice)
        self.no_search_labels()
        self.hit_outputs(hits, qshown, dbshown)

    def t_allpairs_global(self):
        a = self.a
        src = self.filtered(self.inp, 1)
        mode = a.opts.get('qmask', 'dust')
        cmds = self.masking('qmask', 'standard')
        f = self.masked(src, cmds)
        shown = src if mode == 'none' else f
        args = []
        if a.flag('acceptall'):
            args.append('--accept-all')
        for v in ('maxaccepts', 'maxrejects', 'wordlength', 'minwordmatches'):
            a.take(v)  # vsearch's do nothing here (decision 29)
        args += self.search_args(words=False)
        h = self.step('allpr', f, args, name='uc')
        self.no_search_labels()
        self.hit_outputs(lambda kind: h, shown, shown)

    def t_orient(self):
        a = self.a
        qmode = a.take('qmask', 'dust')
        if qmode not in ('dust', 'soft', 'none'):
            raise TranslationError(f'--qmask {qmode}')
        q = self.inp
        if qmode == 'none':
            for o in ('fastaout', 'fastqout', 'notmatched'):
                if a.has(o):
                    raise TranslationError(
                        f'--{o} with --qmask none: lcmask --clear upper-cases '
                        'the queries vsearch writes as they are')
            q = self.chain(q, [('lcmask', ['--clear'])])
        db = a.take('db')
        if db is None:
            raise TranslationError('--orient needs --db')
        if sniff(db) == 'udb':
            self.lengths(1)
            a.take('dbmask')
            a.take('hardmask')
        else:
            mode = a.take('dbmask', 'dust')
            hard = a.flag('hardmask')
            cmds = {('dust', False): [('lcmask', [])],
                    ('dust', True): [('lcmask', ['--hard'])],
                    ('soft', False): [],
                    ('soft', True): [('lcmask', ['--hard', '--no-dust'])],
                    ('none', False): [('lcmask', ['--clear'])],
                    ('none', True): [('lcmask', ['--clear'])]}.get((mode, hard))
            if cmds is None:
                raise TranslationError(f'--dbmask {mode}')
            args = self.length_filter(1, what='database')
            s = db
            if args is not None:
                s = self.run('sqfilt', ['-i', db] + args, name='db')
            db = self.chain(s, cmds, first_reads_file=(s == db)) if cmds else s
        args = ['--db', db]
        if a.has('wordlength'):
            args += ['--word-len', str(a.take('wordlength'))]
        notm = self.tmp('notmatched')
        args += ['--not-oriented', notm]
        if a.has('tabbedout'):
            args += ['--tabbedout', self.path('tabbedout')]
        if a.has('log'):
            args += ['--report', self.path('log')]
        out = self.step('dbornt', q, args, name='oriented')
        if a.has('tabbedout'):
            self.output('tabbedout', self.path('tabbedout'))
        if a.has('log'):
            self.output('log', self.path('log'), 'orient-log')
        self.no_search_labels()
        fmt = sniff(self.inp)
        self.records(out, fmt, 'fastaout', 'fastqout')
        if a.has('notmatched'):
            if fmt == 'fastq':
                self.output('notmatched', notm)
            else:
                self.fasta('notmatched', notm, fmt)

    # clustering

    def cluster(self, sort: Optional[str], unoise: bool = False):
        a = self.a
        extra = []
        if unoise:
            extra = ['--min-size', str(a.take('minsize', 8))]
            if extra[1] == '1':
                extra = []
        src = self.filtered(self.inp, 32, extra)
        if extra and src == self.inp:
            raise TranslationError('--minsize with --minseqlength 0')
        cmds = self.masking('qmask', 'standard')
        if sort:
            cmds = cmds + [('sortby', [sort])]
        s = self.masked(src, cmds)
        args = []
        if unoise:
            a.take('id')  # vsearch ignores --id for --cluster_unoise
            if a.has('weak_id'):
                args += ['--id', str(a.take('weak_id'))]
            if a.has('unoise_alpha'):
                args += ['--alpha', str(a.take('unoise_alpha'))]
        elif not a.has('id'):
            raise TranslationError(f'--{self.cmd} needs --id')
        if self.cmd == 'cluster_fast' and not a.has('maxrejects'):
            args += ['--maxrejects', '8']
        if self.cmd == 'cluster_smallmem' and not a.flag('usersort'):
            self.caveat('without --usersort vsearch requires the input sorted '
                        'by length')
        args += self.search_args(hitchoice=False, skip=('weak_id', 'id') if unoise else ())
        if a.flag('sizeorder'):
            args.append('--size-order')
        sizein, sizeout = a.flag('sizein'), a.flag('sizeout')
        sz = ([] if sizein else ['--no-sizein']) + \
            ([] if sizeout else ['--no-sizeout'])
        hits = self.tmp('uc')
        args += sz + ['--hits', hits]
        cent = self.step('unoise' if unoise else 'clust', s, args,
                         name='centroids')
        self.no_labels(f'relabelling the centroids of --{self.cmd} is not '
                       'translated', skip=('sizeout',))
        for o in ('alnout', 'blast6out', 'userout', 'samout', 'fastapairs',
                  'qsegout', 'tsegout', 'matched', 'notmatched', 'otutabout',
                  'mothur_shared_out', 'biomout', 'centroid_sizeout',
                  'clusterout_id', 'clusterout_sort', 'cons_truncate'):
            if a.has(o):
                raise TranslationError(f'--{o} with --{self.cmd}: no rsearch '
                                       'translation')
        fmt = sniff(self.inp)
        self.fasta('centroids', cent, fmt)
        if a.has('uc'):
            self.output('uc', hits,
                        'uc-fast' if self.cmd == 'cluster_fast' else 'bytes')
        if a.has('clusters'):
            if self.width != 0 and fmt != 'fastq':
                raise TranslationError(
                    f'--clusters with --fasta_width {self.width}: clspl writes '
                    'records unwrapped (give --fasta_width 0)')
            if sizeout:
                self.caveat('--clusters with --sizeout: vsearch rewrites the '
                            'abundances clspl writes as they are (decision 33)')
            prefix = self.path('clusters')
            self.run('clspl', ['--hits', hits, '-i', s, '--prefix', prefix],
                     name='clspl')
            self.output('clusters', prefix, 'cluster-files')
        for opt, flag in (('msaout', []), ('consout', ['--consensus']),
                          ('profile', ['--profile'])):
            if not a.has(opt):
                continue
            if opt != 'consout':
                self.unwrapped_only(opt)
            out = self.run('clmsa', ['--hits', hits, '-i', s] + flag + sz,
                           name=opt)
            if opt == 'consout':
                self.fasta(opt, out, 'fasta')
            else:
                self.output(opt, out)

    def t_cluster_fast(self):
        self.cluster('--length')

    def t_cluster_size(self):
        self.cluster('--size')

    def t_cluster_smallmem(self):
        self.cluster(None)

    def t_cluster_unoise(self):
        self.cluster('--size', unoise=True)

    # chimeras

    def uchime(self, denovo: bool, uchime2: bool = False, abskew=None):
        a = self.a
        if denovo:
            src = self.filtered(self.inp, 1)
            cmds = self.masking('qmask', 'standard') + [('sortby', ['--size'])]
            s = self.masked(src, cmds)
            args = []
        else:
            s, _ = self.queries()
            db, _ = self.database(1, 'standard')
            args = ['--db', db]
            a.take('abskew')  # used de novo only
            if a.has('strand') and a.opts['strand'] == 'plus':
                a.take('strand')
        if uchime2:
            args.append('--uchime2')
            for v in ('minh', 'mindiv', 'mindiffs'):
                a.take(v)  # vsearch ignores them for UCHIME2
        if denovo:
            skew = a.take('abskew', abskew)
            if skew is not None:
                args += ['--abskew', str(skew)]
        for v in ('minh', 'mindiv', 'mindiffs', 'xn', 'dn'):
            if a.has(v):
                args += ['--' + v, str(a.take(v))]
        a.take('sizein')
        self.no_labels('uchime writes its records unchanged (decision 35)')
        outs = {}
        for v, r in (('chimeras', '--chimeras'), ('borderline', '--borderline')):
            if a.has(v):
                outs[v] = self.tmp(v)
                args += [r, outs[v]]
        for v, r in (('uchimeout', '--table'), ('uchimealns', '--alignments')):
            if a.has(v):
                args += [r, self.path(v)]
        if a.has('alignwidth'):
            args += ['--width', str(a.take('alignwidth'))]
        non = self.step('uchime', s, args, name='nonchimeras')
        for v in ('uchimeout', 'uchimealns'):
            if a.has(v):
                self.output(v, self.path(v))
        fmt = sniff(self.inp)
        self.fasta('nonchimeras', non, fmt)
        for v, f in outs.items():
            self.fasta(v, f, fmt)

    def t_uchime_denovo(self):
        self.uchime(True)

    def t_uchime2_denovo(self):
        self.uchime(True, uchime2=True)

    def t_uchime3_denovo(self):
        self.uchime(True, uchime2=True, abskew=16)

    def t_uchime_ref(self):
        self.uchime(False)

    def t_chimeras_denovo(self):
        a = self.a
        src = self.filtered(self.inp, 1)
        cmds = self.masking('qmask', 'standard') + [('sortby', ['--size'])]
        s = self.masked(src, cmds)
        args = []
        for v, r in (('chimeras_parts', '--parts'),
                     ('chimeras_parents_max', '--max-parents'),
                     ('chimeras_length_min', '--min-segment'),
                     ('chimeras_diff_pct', '--diff-pct'),
                     ('abskew', '--abskew'), ('alignwidth', '--width')):
            if a.has(v):
                args += [r, str(a.take(v))]
        a.take('sizein')
        self.no_labels('mchime writes its records unchanged (decision 36)')
        chim = None
        if a.has('chimeras'):
            chim = self.tmp('chimeras')
            args += ['--chimeras', chim]
        for v, r in (('tabbedout', '--table'), ('alnout', '--alignments')):
            if a.has(v):
                args += [r, self.path(v)]
        non = self.step('mchime', s, args, name='nonchimeras')
        for v in ('tabbedout', 'alnout'):
            if a.has(v):
                self.output(v, self.path(v))
        fmt = sniff(self.inp)
        self.fasta('nonchimeras', non, fmt)
        if chim:
            self.fasta('chimeras', chim, fmt)

    # taxonomy

    def t_sintax(self):
        a = self.a
        db = a.take('db')
        if db is None:
            raise TranslationError('--sintax needs --db')
        if sniff(db) != 'udb':
            mode = a.take('dbmask', 'dust')
            if mode not in ('dust', 'soft', 'none'):
                raise TranslationError(f'--dbmask {mode}')
            args = self.length_filter(32, what='database')
            s = db
            if args is not None:
                s = self.run('sqfilt', ['-i', db] + args, name='db')
            if mode == 'none':
                s = self.chain(s, [('lcmask', ['--clear'])],
                               first_reads_file=(s == db))
            db = s
        else:
            self.lengths(32)
            a.take('dbmask')
        q = self.chain(self.inp, [('lcmask', ['--clear'])])
        args = ['--db', db]
        for v, r in (('wordlength', '--word-len'), ('strand', '--strand'),
                     ('sintax_cutoff', '--cutoff')):
            if a.has(v):
                args += [r, str(a.take(v))]
        if a.flag('sintax_random'):
            args.append('--random')
        args += self.seed()
        if a.has('log'):
            args += ['--report', self.path('log')]
        out = self.run('sintax', args, stdin=q, name='sintax')
        self.fastq('tabbedout', out)
        if a.has('log'):
            self.output('log', self.path('log'), 'sintax-log')

    # UDB files

    def t_makeudb_usearch(self):
        a = self.a
        src = self.filtered(self.inp, 32)
        s = self.masked(src, self.masking('dbmask', 'udb'))
        args = []
        if a.has('wordlength'):
            args += ['--word-len', str(a.take('wordlength'))]
        out = self.step('mkudb', s, args, name='udb')
        self.fastq('output', out)

    def t_udbinfo(self):
        self.log_output('udbinfo', ['-i', self.inp], 'udbinfo-log')

    def t_udbstats(self):
        self.log_output('udbstat', ['-i', self.inp], 'udbstats-log')


def translate(argv: List[str], workdir: str = '.') -> Pipeline:
    """The rsearch pipeline for vsearch's command line `argv` (with or
    without the program name), writing its files into `workdir`."""
    return _Translator(list(argv), workdir).translate()
