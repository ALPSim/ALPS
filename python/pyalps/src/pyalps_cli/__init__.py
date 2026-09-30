# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Shell entry points for bundled ALPS programs, without importing pyalps.

Keep this outside pyalps: importing that package also loads the scientific
Python stack, which the native command-line applications do not need.
"""

import os
from pathlib import Path
import sys


def _run(program):
    package = Path(__file__).resolve().parent.parent / "pyalps"
    executable = package / "bin" / program
    if not executable.is_file():
        print(
            f"{program}: this pyalps installation does not bundle the executable. "
            "Install a wheel built with PYALPS_BUNDLE_APPLICATIONS=ON, or invoke "
            "the executable from your ALPS SDK by its full path.",
            file=sys.stderr,
        )
        return 127

    env = os.environ.copy()
    # Respect explicit resource overrides, just like pyalps.tools.
    env.setdefault("ALPS_XML_PATH", str(package / "xml"))
    env.setdefault("ALPS_BIN_PATH", str(package / "bin"))
    try:
        # Use the bundled absolute path, never PATH: it may contain this
        # very launcher or an unrelated ALPS install. Replace the process so
        # arguments, streams, exit status and signals reach the native tool.
        os.execve(str(executable), [str(executable), *sys.argv[1:]], env)
    except OSError as error:
        print(f"{program}: cannot execute {executable}: {error}", file=sys.stderr)
        return 126


def checksign():
    return _run("checksign")


def dirloop_sse():
    return _run("dirloop_sse")


def dmft():
    return _run("dmft")


def dmrg():
    return _run("dmrg")


def fulldiag():
    return _run("fulldiag")


def fulldiag_evaluate():
    return _run("fulldiag_evaluate")


def hirschfye():
    return _run("hirschfye")


def hybridization():
    return _run("hybridization")


def interaction():
    return _run("interaction")


def loop():
    return _run("loop")


def qwl():
    return _run("qwl")


def qwl_evaluate():
    return _run("qwl_evaluate")


def simplemc():
    return _run("simplemc")


def sparsediag():
    return _run("sparsediag")


def spinmc():
    return _run("spinmc")


def spinmc_evaluate():
    return _run("spinmc_evaluate")


def worm():
    return _run("worm")


def worm_evaluate():
    return _run("worm_evaluate")


def parameter2xml():
    return _run("parameter2xml")


def printgraph():
    return _run("printgraph")
