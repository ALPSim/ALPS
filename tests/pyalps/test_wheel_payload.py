#!/usr/bin/env python3
# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Guard what the wheel actually ships next to the extension modules.

The wheel bundles the ALPS application programs (pyalps/bin) and the ALPS
shared libraries they link against (pyalps/lib), while auditwheel/delocate
vendor the genuinely external dependencies alongside them (pyalps.libs on
Linux, pyalps/.dylibs on macOS).  Neither the ordinary binding tests nor the
wheel smoke tests ever load one of those programs, so two mistakes used to
travel undetected: shipping the same library more than once, and shipping
programs whose dynamic dependencies cannot be resolved from inside the
installed package.  These tests run after `repair-wheel-command`, so they
check the repaired artifact rather than the build tree.
"""

from __future__ import annotations

import collections
from pathlib import Path
import re
import subprocess

import pytest

import pyalps


EXPECTED_BUNDLED_PROGRAMS = {
    "checksign",
    "dirloop_sse",
    "dmft",
    "dmrg",
    "fulldiag",
    "fulldiag_evaluate",
    "hirschfye",
    "hybridization",
    "interaction",
    "loop",
    "qwl",
    "qwl_evaluate",
    "simplemc",
    "sparsediag",
    "spinmc",
    "spinmc_evaluate",
    "worm",
    "worm_evaluate",
}


def _package_dir() -> Path:
    return Path(pyalps.__file__).resolve().parent


def _library_dirs() -> list[Path]:
    """Every directory in the installed package that holds bundled libraries."""
    pkg = _package_dir()
    candidates = [pkg / "lib", pkg / ".dylibs", pkg.parent / f"{pkg.name}.libs"]
    return [d for d in candidates if d.is_dir()]


def _library_stem(name: str) -> str:
    """Reduce a shared-library file name to the library it is a copy of.

    ``libalps.so.2.3.4``, ``libalps-034f2e8c.so.2.3.4`` and ``libalps.dylib``
    all reduce to ``libalps``: the version suffixes and the content hash that
    auditwheel/delocate append are what make duplicate copies look distinct.
    """
    if ".so" in name:
        stem = name.split(".so", 1)[0]
    elif name.endswith(".dylib"):
        stem = re.sub(r"(\.\d+)+$", "", name[: -len(".dylib")])
    else:
        return ""
    return re.sub(r"-[0-9a-f]{6,}$", "", stem)


def test_no_shared_library_is_bundled_twice():
    library_dirs = _library_dirs()
    if not library_dirs:
        pytest.skip("no bundled libraries in this install (source tree or SDK install)")

    seen: dict[str, list[Path]] = collections.defaultdict(list)
    for directory in library_dirs:
        for entry in sorted(directory.iterdir()):
            if not entry.is_file():
                continue
            stem = _library_stem(entry.name)
            if stem:
                seen[stem].append(entry)

    pkg_parent = _package_dir().parent
    duplicates = {
        stem: [str(p.relative_to(pkg_parent)) for p in paths]
        for stem, paths in seen.items()
        if len(paths) > 1
    }
    assert not duplicates, (
        "the same shared library is bundled more than once; each copy is dead "
        f"weight in the wheel: {duplicates}"
    )


def test_every_bundled_program_can_be_loaded():
    bin_dir = _package_dir() / "bin"
    if not bin_dir.is_dir():
        pytest.skip("this install does not bundle the ALPS programs")

    programs = sorted(p for p in bin_dir.iterdir() if p.is_file())
    assert programs, f"{bin_dir} exists but is empty"
    assert {p.name for p in programs} == EXPECTED_BUNDLED_PROGRAMS

    # Signatures the dynamic loader emits when a dependency cannot be resolved
    # from inside the installed package.  A program is free to reject --help
    # however it likes -- several of these tools have no option parsing and
    # abort on an uncaught C++ exception, so the exit status alone says nothing
    # -- but it is not free to fail to start.  dyld and glibc/musl both print a
    # distinctive message before dying, which is what this matches on.
    loader_errors = (
        "error while loading shared libraries",  # glibc
        "cannot open shared object file",        # glibc, detail line
        "Error loading shared library",          # musl
        "Library not loaded",                    # dyld
        "image not found",                       # dyld
        "ymbol not found",                       # dyld, either capitalisation
    )

    failures = []
    for program in programs:
        try:
            proc = subprocess.run(
                [str(program), "--help"],
                capture_output=True,
                text=True,
                timeout=60,
            )
        except subprocess.TimeoutExpired:
            # It started running, which is the only thing under test here.
            continue
        except OSError as exc:
            failures.append(f"{program.name}: {exc}")
            continue
        output = f"{proc.stdout}\n{proc.stderr}"
        hit = next((sig for sig in loader_errors if sig in output), None)
        if hit is not None:
            failures.append(f"{program.name}: loader error ({hit!r})")
        elif proc.returncode == 127:
            failures.append(
                f"{program.name}: exited 127: {output.strip()[:200]}"
            )

    assert not failures, "bundled programs that cannot start:\n  " + "\n  ".join(failures)


def test_version_is_inherited_from_the_repository():
    """pyalps' version must come from ALPS_VERSION.txt, not a second copy.

    The numeric core is read from that file at build time; a prerelease label,
    which cannot live there because project(VERSION ...) rejects a non-numeric
    version, comes from the ALPS_VERSION_PRERELEASE environment variable. So the
    installed version must be the file's contents followed by nothing or by a
    PEP 440 prerelease segment -- never a different number.
    """
    import pyalps

    version_file = Path(__file__).resolve().parents[2] / "ALPS_VERSION.txt"
    if not version_file.is_file():
        pytest.skip("not running from an ALPS checkout")

    core = version_file.read_text(encoding="utf-8").strip().splitlines()[0].strip()
    assert re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", core), core

    installed = pyalps.__version__
    assert installed.startswith(core), (
        f"pyalps.__version__ is {installed!r} but ALPS_VERSION.txt says {core!r}; "
        "the version is no longer inherited from the repository"
    )
    suffix = installed[len(core):]
    assert re.fullmatch(r"|(a|b|rc)[0-9]+|\.dev[0-9]+", suffix), (
        f"unexpected version suffix {suffix!r}"
    )
