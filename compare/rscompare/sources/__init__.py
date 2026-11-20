"""Case sources. Each source writes its inputs with a generator of
`compare/gen` from a seed, exactly as the zsh runner did, and returns
its cases: vsearch command lines whose outputs go under `OUT`.

`SOURCES[name]` is a `Source`: the generator, the default seeds (those
the zsh runner passed, from compare/README.md) and the function making
the cases.
"""

import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from typing import Callable, Dict, List

from ..cases import OUT, Case

GEN = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), 'gen')


@dataclass
class Source:
    name: str
    seeds: List[int]
    make: Callable[[str, int, 'Env'], List[Case]]
    about: str = ''


@dataclass
class Env:
    bindir: str
    vsearch: str


SOURCES: Dict[str, Source] = {}


def source(name: str, seeds, about: str = ''):
    def wrap(fn):
        SOURCES[name] = Source(name, list(seeds), fn, about)
        return fn
    return wrap


def gen(script: str, d: str, seed, *extra):
    subprocess.run([sys.executable, os.path.join(GEN, script + '.py'), d,
                    str(seed)] + [str(x) for x in extra], check=True,
                   stdout=subprocess.DEVNULL)


def argv(text: str, **kw) -> List[str]:
    """A command line from text, with {O} the output directory."""
    return shlex.split(text.format(O=OUT, **kw))


def read(path: str) -> str:
    with open(path) as f:
        return f.read()


def seeds(a: int, b: int) -> List[int]:
    return list(range(a, b + 1))


from . import preprocess, analysis, fuzz, datasets  # noqa: E402,F401
