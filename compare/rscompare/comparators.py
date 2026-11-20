"""Named comparators: how a vsearch output, error or warning is compared
with rsearch's. Each normalisation the zsh runners made lives here once.

A file comparator takes the vsearch and rsearch files and returns None
when they agree, or a short description of the first difference. A
missing file reads as empty, as the runners did with `touch`.
"""

import os
import re
from typing import Callable, Dict, List, Optional


def _read(path: str) -> bytes:
    try:
        with open(path, 'rb') as f:
            return f.read()
    except OSError:
        return b''


def _first_diff(v: bytes, r: bytes) -> Optional[str]:
    if v == r:
        return None
    vl, rl = v.split(b'\n'), r.split(b'\n')
    for k, (a, b) in enumerate(zip(vl, rl)):
        if a != b:
            return (f'line {k + 1}: vsearch {a[:120]!r} | rsearch {b[:120]!r}')
    return f'lengths differ: vsearch {len(vl)} lines, rsearch {len(rl)} lines'


def _lines_view(fn: Callable[[List[bytes]], List[bytes]]):
    """A comparator that compares the files after a normalisation of their
    lines (`fn` on each side)."""
    def compare(v: str, r: str) -> Optional[str]:
        return _first_diff(b'\n'.join(fn(_read(v).split(b'\n'))),
                           b'\n'.join(fn(_read(r).split(b'\n'))))
    return compare


def _sides(vfn, rfn=None):
    """A comparator normalising the vsearch file with vfn and the rsearch
    file with rfn (or as it is)."""
    def compare(v: str, r: str) -> Optional[str]:
        a = vfn(_read(v))
        b = rfn(_read(r)) if rfn else _read(r)
        return _first_diff(a, b)
    return compare


def bytes_(v: str, r: str) -> Optional[str]:
    return _first_diff(_read(v), _read(r))


_PERFECT_BUT_ENDS = re.compile(rb'(\d*[ID])?\d*M(\d*[ID])?')


def uc_fast(v: str, r: str) -> Optional[str]:
    """--cluster_fast writes = for an H record's alignment that is perfect
    but for terminal gaps, which clust writes whole (decision 31)."""
    vl, rl = _read(v).split(b'\n'), _read(r).split(b'\n')
    if len(vl) != len(rl):
        return f'lengths differ: vsearch {len(vl)} lines, rsearch {len(rl)}'
    for k, (a, b) in enumerate(zip(vl, rl)):
        x, y = a.split(b'\t'), b.split(b'\t')
        if len(x) > 7 and len(y) > 7 and x[0] == b'H' and x[7] == b'=' \
                and _PERFECT_BUT_ENDS.fullmatch(y[7]):
            y[7] = b'='
        if x != y:
            return f'line {k + 1}: vsearch {a[:120]!r} | rsearch {b[:120]!r}'
    return None


def _biom(lines):
    """BIOM without the lines naming the file, the program and the date
    (decision 26)."""
    return [l for l in lines if not (l.startswith(b'\t"id":"')
                                     or b'"generated_by"' in l
                                     or b'"date"' in l)]


def log_report(data: bytes) -> bytes:
    """The report in a vsearch --log file, without the lines the log adds
    before (program, command line, start time) and after (finish time,
    elapsed time, memory); compare/gen/vsearch_log_report.py."""
    if not data:
        return b''
    lines = data.split(b'\n')[3:]
    ends = [k for k, l in enumerate(lines) if l.startswith(b'Finished')]
    if not ends:
        return b''
    return b'\n'.join(lines[:ends[-1] - 1]) + b'\n'


def _range(start: bytes, end: bytes, anchored: bool):
    """sed -n '/start/,/end/p': from a line matching start to the next line
    matching end (anchored: the line starts with the pattern)."""
    hit = (lambda l, p: l.startswith(p)) if anchored else (lambda l, p: p in l)

    def view(data: bytes) -> bytes:
        out, on = [], False
        for l in data.split(b'\n'):
            if not on:
                if hit(l, start):
                    on = True
                    out.append(l + b'\n')
                continue
            out.append(l + b'\n')
            if hit(l, end):
                on = False
        return b''.join(out)
    return view


def _grep(prefix: bytes):
    def view(data: bytes) -> bytes:
        return b''.join(l + b'\n' for l in data.split(b'\n')
                        if l.startswith(prefix))
    return view


def udbstats_log(data: bytes) -> bytes:
    """--udbstats' report in --log: from the database summary (line 4) to
    the line before the blank line and Finished."""
    if not data:
        return b''
    lines = data.split(b'\n')
    if lines and lines[-1] == b'':
        lines = lines[:-1]
    keep = lines[3:len(lines) - 4]
    return b''.join(l + b'\n' for l in keep)


def _alnout(data: bytes) -> bytes:
    """--alnout without its two opening lines (decision 25)."""
    return b'\n'.join(data.split(b'\n')[2:]) if data else b''


_UR = re.compile(rb'\tUR:file:\S+')


def _sam(data: bytes) -> bytes:
    """SAM without @PG lines, and the database file in @SQ's UR field."""
    return b'\n'.join(_UR.sub(b'\tUR:file:', l) for l in data.split(b'\n')
                      if not l.startswith(b'@PG'))


def upper_seq(data: bytes) -> bytes:
    """Sequence lines upper-cased, FASTA or FASTQ (the input of
    lcmask --clear's comparison with --qmask none)."""
    lines = data.split(b'\n')
    fq = data[:1] == b'@'
    return b'\n'.join(l.upper() if (k % 4 == 1 if fq else not l.startswith(b'>'))
                      else l for k, l in enumerate(lines))


def fastq_to_fasta(width: int):
    """FASTQ as vsearch writes it in FASTA at a width (0: unwrapped)."""
    def view(data: bytes) -> bytes:
        lines = data.split(b'\n')
        out = []
        for k in range(0, len(lines) - 3, 4):
            out.append(b'>' + lines[k][1:] + b'\n')
            s = lines[k + 1]
            if width:  # an empty sequence: no line, as vsearch wraps
                out += [s[i:i + width] + b'\n' for i in range(0, len(s), width)]
            else:
                out.append(s + b'\n')
        return b''.join(out)
    return view


def cluster_files(v: str, r: str) -> Optional[str]:
    """The files of --clusters PREFIX: the same names after the prefix,
    each byte for byte."""
    def files(prefix):
        d, base = os.path.split(prefix)
        try:
            names = os.listdir(d or '.')
        except OSError:
            return {}
        return {n[len(base):]: os.path.join(d, n) for n in names
                if n.startswith(base)}
    a, b = files(v), files(r)
    if set(a) != set(b):
        return (f'files differ: vsearch {len(a)} ({sorted(a)[:5]}), '
                f'rsearch {len(b)} ({sorted(b)[:5]})')
    for k in sorted(a, key=lambda x: (len(x), x)):
        d = bytes_(a[k], b[k])
        if d:
            return f'cluster {k}: {d}'
    return None


FILES: Dict[str, Callable[[str, str], Optional[str]]] = {
    'bytes': bytes_,
    'uc-fast': uc_fast,
    'biom': _lines_view(_biom),
    'log-report': _sides(log_report),
    'orient-log': _sides(_range(b'Forward oriented', b'Total number', True)),
    'sintax-log': _sides(_grep(b'Classified')),
    'udbinfo-log': _sides(_range(b'Seqs', b'DBAccel', False)),
    'udbstats-log': _sides(udbstats_log),
    'alnout': _sides(_alnout),
    'sam': _sides(_sam, _sam),
    'upper-seq': _sides(upper_seq),
    'fastq-to-fasta': _sides(fastq_to_fasta(80)),
    'fastq-to-fasta0': _sides(fastq_to_fasta(0)),
    'cluster-files': cluster_files,
}


# error messages: vsearch's first "Fatal error: X" against rsearch's first
# "Error: X", each normalised by a kind's substitutions

_V_FATAL = re.compile(r'Fatal error: (.*)')
_R_ERROR = re.compile(r'^Error: (.*)', re.M)


def vsearch_error(stderr: str) -> str:
    m = _V_FATAL.search(stderr)
    return m.group(1) if m else ''


def rsearch_error(stderr: str) -> str:
    m = _R_ERROR.search(stderr)
    return m.group(1) if m else ''


_ABUND = [(r'.*abundance.*', 'abundance')]
_UDB_V = [(r'Unable to read.*', 'truncated'), (r'Invalid UDB file', 'invalid'),
          (r'Incorrect UDB file size', 'trailing'), (r'.*too long.*', 'too long')]
_UDB_R = [(r'.*(shorter than its header|truncated)', 'truncated'),
          (r'.*data after its end', 'trailing'), (r'^invalid UDB file$', 'invalid'),
          (r'.*too long.*', 'too long')]

# kind -> (vsearch substitutions, rsearch substitutions)
ERRORS: Dict[str, tuple] = {
    'exact': ([], []),
    'abundance': (_ABUND, _ABUND),
    'joinpa': ([], [(r'^invalid \(([a-z]+)\) abundance annotation in read .*',
                     r'Invalid (\1) abundance annotation in FASTA file header')]),
    'mergpa': (_ABUND + [(r'FASTQ quality value .*(above|below).*', r'\1')],
               _ABUND + [(r'quality score .*(above|below).*', r'\1')]),
    'mkudb': (_ABUND + [(r'.*(empty database|no sequences).*', 'empty')],
              _ABUND + [(r'.*(empty database|no sequences).*', 'empty')]),
    'mkudb-exact': ([], [(r'no sequences: a UDB file needs at least one',
                          'Cannot write a UDB file for an empty database')]),
    'udb': (_UDB_V, _UDB_R),
}


def _apply(msg: str, subs) -> str:
    for pat, rep in subs:
        msg = re.sub(pat, rep, msg)
    return msg


def error_kinds(kind: str, vstderr: str, rstderr: str):
    """The two messages as the kind compares them."""
    v, r = vsearch_error(vstderr), rsearch_error(rstderr)
    if kind == 'invalid':  # udbinfo: both say "invalid", or neither does
        f = lambda m: 'invalid' if 'invalid' in m.lower() else ''
        return f(v), f(r)
    vs, rs = ERRORS[kind]
    return _apply(v, vs), _apply(r, rs)


# warnings: a value extracted from each side's stderr, compared

def _first(pat: str):
    rx = re.compile(pat)

    def get(text: str) -> str:
        m = rx.search(text)
        return m.group(0) if m else ''
    return get


def _count(pat: str):
    return lambda text: str(len(re.findall(pat, text)))


def _nowords(v_phrase: str, r_phrase: str):
    def pair(vtext, rtext):
        def get(text, phrase):
            m = re.search(r'(\d+ of \d+) sequences ' + phrase, text)
            return m.group(1) if m else ''
        return get(vtext, v_phrase), get(rtext, r_phrase)
    return pair


WARNINGS: Dict[str, Callable[[str, str], tuple]] = {
    'stripped': lambda v, r: (_first(r'invalid characters stripped.*')(v),
                              _first(r'invalid characters stripped.*')(r)),
    'missing-abundance': lambda v, r: (_count('Missing abundance')(v),
                                       _count('Missing abundance')(r)),
    'nowords': _nowords('yielded no k-mer', 'give no word'),
    'wordlength-adjusted': lambda v, r: (_count('Wordlength adjusted')(v),
                                         _count('word length of the UDB file')(r)),
}
