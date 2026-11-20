"""Runs vsearch and rsearch pipelines."""

import os
import subprocess
import time
from dataclasses import dataclass, field
from typing import List, Optional

from .translate import Step

TIMEOUT = 1800


@dataclass
class Run:
    exit: int
    stderr: str                  # the failing step's (or vsearch's) stderr
    all_stderr: str              # every step's stderr
    seconds: float
    failed: Optional[str] = None  # the first failing step, rendered
    steps: List[dict] = field(default_factory=list)


def _proc(argv, stdin_path, stdout_path, cwd):
    fin = open(stdin_path, 'rb') if stdin_path else subprocess.DEVNULL
    fout = open(stdout_path, 'wb') if stdout_path else subprocess.DEVNULL
    try:
        p = subprocess.run(argv, stdin=fin, stdout=fout, stderr=subprocess.PIPE,
                           cwd=cwd, timeout=TIMEOUT)
        code, err = p.returncode, p.stderr
    except subprocess.TimeoutExpired as e:
        code, err = -9, (e.stderr or b'') + b'\nrscompare: timed out\n'
    finally:
        if stdin_path:
            fin.close()
        if stdout_path:
            fout.close()
    if code < 0:
        code = 128 - code  # killed by a signal, as a shell reports it
    return code, err.decode('utf-8', 'replace')


def run_vsearch(vsearch: str, argv: List[str], cwd: str) -> Run:
    t = time.time()
    code, err = _proc([vsearch] + argv, None, None, cwd)
    return Run(code, err, err, time.time() - t)


def run_steps(steps: List[Step], bindir: str, cwd: str, outdir: str) -> Run:
    """Runs the steps in order. A failed step's output, and whatever is made
    from it, is incomplete: a later step still reads it on its standard
    input, as a shell pipeline would, but a step that takes it as a file
    argument (a database, a hit file, a second input) is skipped. The exit
    code is the first failing step's."""
    t0 = time.time()
    first = None
    errs, info = [], []
    made, tainted = set(), set()
    prefix = os.path.join(outdir, '')
    for s in steps:
        args_in = [a for a in s.args if a in made]
        no_stdin = s.stdin and not os.path.exists(
            s.stdin if os.path.isabs(s.stdin) else os.path.join(cwd, s.stdin))
        if no_stdin or any(a in tainted for a in args_in):
            info.append({'tool': s.tool, 'skipped': True})
            skipped = [a for a in s.args if a.startswith(prefix)
                       and a not in made] + ([s.stdout] if s.stdout else [])
            made.update(skipped)
            tainted.update(skipped)
            continue
        outs = [a for a in s.args if a.startswith(prefix) and a not in made]
        if s.stdout:
            outs.append(s.stdout)
        t = time.time()
        code, err = _proc(s.argv(bindir), s.stdin, s.stdout, cwd)
        info.append({'tool': s.tool, 'exit': code,
                     'seconds': round(time.time() - t, 4)})
        errs.append(err)
        made.update(outs)
        if code != 0 or (s.stdin in tainted):
            tainted.update(outs)
        if code != 0 and first is None:
            first = (code, err, s)
    total = time.time() - t0
    if first is None:
        return Run(0, '', ''.join(errs), total, None, info)
    code, err, s = first
    return Run(code, err, ''.join(errs), total,
               ' '.join([s.tool] + s.args), info)
