"""Real datasets: the manifest compare/datasets.json, downloaded by
compare/fetch_datasets.py into a data directory ($RSCOMPARE_DATA, by
default $TMPDIR/rsearch-data, one directory per dataset).

The source reads each manifest entry's vsearch command lines and runs
them as cases; the seed is the entry's place in the manifest (from 1).
An entry is a JSON object; the source uses:

- its name: "id", else "name", else "accession";
- its command lines: "vsearch" or "commands" (or "vsearch_commands"), a
  list of strings or of argv lists, without the program name. In them
  `{dir}` is the dataset's directory (relative file names are taken
  from it too); `{O}` is the output directory. Output options without
  `{O}` are moved there, named by the option.

Every command line gets `--threads 1` and `--quiet`, and `--notrunclabels`
where vsearch would cut labels. An entry whose files are missing is
skipped, naming them.

The manifest itself is written by another session; nothing here writes
it.
"""

import json
import os
import shlex

from ..cases import OUT, Case
from ..translate import WHOLE_HEADERS
from ..vsearch_options import Args, TranslationError
from . import GEN, source

MANIFEST = os.path.join(os.path.dirname(GEN), 'datasets.json')

OUTPUTS = set("""
alnout biomout blast6out borderline centroids chimeras clusters consout
dbmatched dbnotmatched eetabbedout fastaout fastaout_discarded
fastaout_discarded_rev fastaout_notmerged_fwd fastaout_notmerged_rev
fastaout_orphans fastaout_orphans_rev fastaout_rev fastapairs fastqout
fastqout_discarded fastqout_discarded_rev fastqout_notmerged_fwd
fastqout_notmerged_rev fastqout_orphans fastqout_orphans_rev fastqout_rev
lcaout log matched mothur_shared_out msaout nonchimeras notmatched
notmatchedfq otutabout output profile qsegout samout tabbedout tsegout uc
uchimealns uchimeout userout
""".split())

# options whose value names an input file
INPUTS = {'db', 'reverse', 'labels', 'label_words'}


def entries():
    if not os.path.exists(MANIFEST):
        return []
    with open(MANIFEST) as f:
        data = json.load(f)
    if isinstance(data, dict):
        data = data.get('datasets', list(data.values()))
    return [e for e in data if isinstance(e, dict)]


def data_dir():
    return os.environ.get('RSCOMPARE_DATA', os.path.join(
        os.environ.get('TMPDIR', '/tmp'), 'rsearch-data'))


def entry_name(e):
    return str(e.get('id') or e.get('name') or e.get('accession') or '?')


def prepare(cmd, d):
    """A manifest command line as a case's argv: inputs from the dataset's
    directory, outputs moved under OUT."""
    words = shlex.split(cmd.replace('{dir}', d)) if isinstance(cmd, str) \
        else [str(w).replace('{dir}', d) for w in cmd]
    if words and not words[0].startswith('-'):
        words = words[1:]  # the program
    a = Args(words)  # checks the options
    out, missing = [], []
    k = 0
    while k < len(words):
        w = words[k]
        name = w.lstrip('-')
        out.append(w)
        if name in OUTPUTS and k + 1 < len(words):
            v = words[k + 1]
            out.append(v.replace('{O}', OUT) if '{O}' in v
                       else os.path.join(OUT, name))
            k += 2
            continue
        if (name == a.command or name in INPUTS) and k + 1 < len(words):
            v = words[k + 1]
            if not os.path.isabs(v):
                v = os.path.join(d, v)
            if not os.path.exists(v):
                missing.append(v)
            out.append(v)
            k += 2
            continue
        k += 1
    have = set(x.lstrip('-') for x in out)
    if 'threads' not in have:
        out += ['--threads', '1']
    if 'quiet' not in have:
        out.append('--quiet')
    if a.command not in WHOLE_HEADERS and 'notrunclabels' not in have:
        out.append('--notrunclabels')
    return out, missing


def _commands(e):
    for key in ('vsearch', 'commands', 'vsearch_commands'):
        if e.get(key):
            v = e[key]
            return v if isinstance(v, list) else [v]
    return []


@source('datasets', [], 'the real datasets of compare/datasets.json, by their '
        'place in it (seeds 1 to N)')
def datasets(T, seed, env):
    es = entries()
    if not 1 <= seed <= len(es):
        return []
    e = es[seed - 1]
    name = entry_name(e)
    d = os.path.join(data_dir(), name)
    cases = []
    for k, cmd in enumerate(_commands(e)):
        try:
            argv, missing = prepare(cmd, d)
        except TranslationError as err:
            argv, missing = [], [f'unparsable command line: {err}']
        if missing:
            why = 'missing ' + ', '.join(missing[:3])
            cases.append(Case(f'{name} #{k}', [], skip=lambda v, r, w=why: w))
            continue
        cases.append(Case(f'{name} #{k}', argv, errors='exact',
                          warnings=['stripped'], coverage={'dataset': name}))
    return cases


def default_seeds():
    return list(range(1, len(entries()) + 1))
