# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Derive the pyalps version from the repository, not from pyproject.toml.

ALPS_VERSION.txt at the repository root is the single source of truth for the
release version; cmake/ALPSVersion.cmake reads the same file to set
ALPS_VERSION_CORE before ``project()``. This provider reads it for the Python
package metadata, so a release bump is one edit rather than two that can drift.

A prerelease label is deliberately *not* in that file: ``project(VERSION ...)``
rejects a non-numeric version, and neither ``find_package()`` matching nor the
library SOVERSION has a notion of prerelease ordering. CMake therefore takes it
from the ALPS_VERSION_PRERELEASE cache variable, set by the release process. The
environment variable of the same name is the equivalent here, using the same
vocabulary, translated to the PEP 440 spelling Python requires:

    (unset)  -> 2.3.4          a final release
    beta.1   -> 2.3.4b1
    alpha.2  -> 2.3.4a2
    rc.1     -> 2.3.4rc1
    dev.3    -> 2.3.4.dev3

Wired up in pyproject.toml as::

    [project]
    dynamic = ["version"]

    [[tool.dynamic-metadata]]
    provider = { path = "_build_support", module = "alps_version" }
"""

from __future__ import annotations

import os
import re
from pathlib import Path

#: Where ALPS_VERSION.txt can be, in the order tried.
#:
#: A build from the repository finds it three levels up. An sdist build finds it
#: at the root, because pyproject.toml force-includes it there -- the sdist is
#: rooted at this project directory and cannot reach outside itself, which is
#: also why the CMake side keeps its own repo-or-_vendor check.
_CANDIDATES = (
    Path("ALPS_VERSION.txt"),
    Path("..") / ".." / ".." / "ALPS_VERSION.txt",
)

_CORE = re.compile(r"^(?P<core>[0-9]+\.[0-9]+\.[0-9]+)$")

#: CMake's prerelease vocabulary mapped to PEP 440 separators.
_PRERELEASE_KINDS = {
    "alpha": "a",
    "a": "a",
    "beta": "b",
    "b": "b",
    "rc": "rc",
    "c": "rc",
    "pre": "rc",
    "preview": "rc",
    "dev": ".dev",
}

_PRERELEASE = re.compile(r"^(?P<kind>[A-Za-z]+)[.\-_]?(?P<number>[0-9]+)?$")


def _read_core() -> str:
    """Return MAJOR.MINOR.PATCH from ALPS_VERSION.txt."""
    for candidate in _CANDIDATES:
        if not candidate.is_file():
            continue
        text = candidate.read_text(encoding="utf-8").strip().splitlines()
        core = text[0].strip() if text else ""
        match = _CORE.match(core)
        if match is None:
            raise RuntimeError(
                f"{candidate} must contain exactly MAJOR.MINOR.PATCH, but reads "
                f"{core!r}. A prerelease label belongs in the "
                f"ALPS_VERSION_PRERELEASE environment variable, and the leading "
                f"'v' of a release tag is not part of the version."
            )
        return match.group("core")

    tried = ", ".join(str(path) for path in _CANDIDATES)
    raise RuntimeError(
        "Cannot find ALPS_VERSION.txt, which supplies the pyalps version. "
        f"Looked in: {tried} (relative to {Path.cwd()}). A build from the ALPS "
        "repository finds it at the repository root; an sdist carries a copy, "
        "placed there by the sdist.force-include entry in pyproject.toml."
    )


def _pep440_suffix(label: str) -> str:
    """Translate a CMake prerelease label into its PEP 440 spelling."""
    label = label.strip()
    if not label:
        return ""

    match = _PRERELEASE.match(label)
    if match is None:
        raise RuntimeError(
            f"ALPS_VERSION_PRERELEASE={label!r} is not a recognised label. "
            f"Use a kind and a number, e.g. beta.1, rc.2 or dev.3; the "
            f"supported kinds are {sorted(set(_PRERELEASE_KINDS))}."
        )

    kind = match.group("kind").lower()
    if kind not in _PRERELEASE_KINDS:
        raise RuntimeError(
            f"ALPS_VERSION_PRERELEASE={label!r} has unknown kind {kind!r}. "
            f"Supported kinds: {sorted(set(_PRERELEASE_KINDS))}."
        )

    return f"{_PRERELEASE_KINDS[kind]}{match.group('number') or '0'}"


def version() -> str:
    """The full PEP 440 version for this build."""
    return _read_core() + _pep440_suffix(os.environ.get("ALPS_VERSION_PRERELEASE", ""))


def dynamic_metadata(settings, project):  # noqa: ARG001 - provider protocol
    """scikit-build-core dynamic-metadata 0.3 hook."""
    if settings:
        raise RuntimeError(
            "The alps_version provider takes no settings; the version comes "
            "from ALPS_VERSION.txt and ALPS_VERSION_PRERELEASE."
        )
    return {"version": version()}


if __name__ == "__main__":  # a convenience for the release process
    print(version())
