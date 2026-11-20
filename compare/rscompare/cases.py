"""A case, and running one: vsearch on its command line, rsearch on the
translated pipeline (or, for a probe, on the steps the case gives), and
every output, the exit code, the error and the warnings compared."""

import os
import shutil
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

from . import comparators as C
from .execute import Run, run_steps, run_vsearch
from .translate import Output, Pipeline, Step, TranslationError, translate

OUT = '@OUT@'   # stands for the case's vsearch output directory in argv


@dataclass
class Case:
    name: str
    vsearch: List[str]                 # argv, outputs under OUT
    errors: str = 'exit'               # 'exit' or a kind of comparators.ERRORS
    warnings: List[str] = field(default_factory=list)
    outputs_on_error: bool = True      # compare outputs when vsearch fails
    exit_map: Dict[Tuple[int, int], int] = field(default_factory=dict)
    # a probe: rsearch steps given by the case instead of the translation,
    # called with the rsearch directory; returns steps and {option: file}
    probe: Optional[Callable[[str], Tuple[List[Step], Dict[str, str]]]] = None
    # comparator per vsearch option, overriding the translation's
    comparators: Dict[str, str] = field(default_factory=dict)
    only: Optional[List[str]] = None   # compare only these outputs
    # outputs and warnings compared only when vsearch succeeds
    success_only: List[str] = field(default_factory=list)
    skip: Optional[Callable[[Run, Run], Optional[str]]] = None
    coverage: Dict[str, object] = field(default_factory=dict)
    # counts taken from vsearch's outputs for the coverage record:
    # {name: (option, kind)} with kind 'records', 'uc' or 'lines'
    counts: Dict[str, Tuple[str, str]] = field(default_factory=dict)


def _subst(argv: List[str], vdir: str) -> List[str]:
    return [a.replace(OUT, vdir) for a in argv]


def _count(path: str, kind: str):
    try:
        data = open(path, 'rb').read()
    except OSError:
        return None
    lines = data.split(b'\n')
    if lines and lines[-1] == b'':
        lines.pop()
    if kind == 'lines':
        return len(lines)
    if kind == 'records':
        if not data:
            return 0
        if data[:1] == b'@':
            return len(lines) // 4
        return sum(1 for l in lines if l.startswith(b'>'))
    if kind == 'uc':
        out = {}
        for l in lines:
            out[l[:1].decode()] = out.get(l[:1].decode(), 0) + 1
        return out
    if kind.startswith('col'):  # col N VALUE: rows whose column N is VALUE
        _, n, value = kind.split(' ', 2)
        return sum(1 for l in lines if len(l.split(b'\t')) > int(n)
                   and l.split(b'\t')[int(n)].decode() == value)
    return None


def run_case(case: Case, workdir: str, vsearch: str, bindir: str,
             keep: bool = False) -> dict:
    if not case.vsearch:  # nothing to run (a dataset not downloaded)
        why = case.skip(None, None) if case.skip else 'no command'
        return {'case': case.name, 'outcome': 'skip', 'reason': why}
    vdir = os.path.join(workdir, 'v')
    rdir = os.path.join(workdir, 'r')
    shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(vdir)
    os.makedirs(rdir)
    argv = _subst(case.vsearch, vdir)
    rec = {'case': case.name, 'vsearch': argv,
           'coverage': dict(case.coverage)}
    try:
        if case.probe:
            steps, outmap = case.probe(rdir)
            pipe = Pipeline('probe', argv, steps,
                            [Output(k, _subst([k], vdir)[0], v)
                             for k, v in outmap.items()])
            for o in pipe.outputs:  # the vsearch file of each option
                o.vsearch = _option_value(argv, o.option)
            rec['kind'] = 'probe'
        else:
            pipe = translate(argv, rdir)
            rec['kind'] = 'translated'
    except TranslationError as e:
        rec.update(outcome='error', reason=f'translation: {e}')
        return rec
    rec['rsearch'] = pipe.render().split('\n')
    rec['caveats'] = pipe.caveats
    vr = run_vsearch(vsearch, argv, workdir)
    rr = run_steps(pipe.steps, bindir, workdir, rdir)
    rec['exit'] = {'vsearch': vr.exit, 'rsearch': rr.exit}
    rec['seconds'] = {'vsearch': round(vr.seconds, 4),
                      'rsearch': round(rr.seconds, 4)}
    if rr.failed:
        rec['failed_step'] = rr.failed
    for name, (opt, kind) in case.counts.items():
        path = _option_value(argv, opt)
        if path:
            rec['coverage'][name] = _count(path, kind)
    verr = C.vsearch_error(vr.stderr)
    rec['coverage']['exit'] = vr.exit
    if verr:
        rec['coverage']['error'] = verr[:60]
    diffs = []
    if case.skip:
        why = case.skip(vr, rr)
        if why:
            rec.update(outcome='skip', reason=why)
            _cleanup(workdir, keep, False)
            return rec
    rexit = case.exit_map.get((vr.exit, rr.exit), rr.exit)
    if vr.exit != rexit:
        diffs.append(f'exit: vsearch {vr.exit}, rsearch {rr.exit}')
    if case.errors != 'exit':
        v, r = C.error_kinds(case.errors, vr.stderr, rr.stderr)
        if v != r:
            diffs.append(f'error: vsearch {v!r} | rsearch {r!r}')
    for w in case.warnings:
        if vr.exit != 0 and w in case.success_only:
            continue
        v, r = C.WARNINGS[w](vr.all_stderr, rr.all_stderr)
        if v != r:
            diffs.append(f'warning {w}: vsearch {v!r} | rsearch {r!r}')
    if vr.exit == 0 or case.outputs_on_error:
        for o in pipe.outputs:
            if case.only is not None and o.option not in case.only:
                continue
            if vr.exit != 0 and o.option in case.success_only:
                continue
            comp = case.comparators.get(o.option, o.comparator)
            d = C.FILES[comp](o.vsearch, o.rsearch)
            if d:
                diffs.append(f'--{o.option} ({comp}): {d}')
    rec['outcome'] = 'fail' if diffs else 'pass'
    if diffs:
        rec['diffs'] = diffs
        rec['stderr'] = {'vsearch': C.vsearch_error(vr.stderr),
                         'rsearch': C.rsearch_error(rr.stderr)}
    _cleanup(workdir, keep, bool(diffs))
    return rec


def _option_value(argv: List[str], option: str) -> Optional[str]:
    for k, a in enumerate(argv[:-1]):
        if a.lstrip('-') == option:
            return argv[k + 1]
    return None


def _cleanup(workdir: str, keep: bool, failed: bool):
    if not keep and not failed:
        shutil.rmtree(workdir, ignore_errors=True)
