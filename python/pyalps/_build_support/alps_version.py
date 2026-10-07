# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Derive the pyalps version and conditional wheel commands.

cmake/ALPS_VERSION.txt in the repository is the single source of truth for the
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
    dynamic = ["version", "scripts"]

    [[tool.dynamic-metadata]]
    provider = { path = "_build_support", module = "alps_version" }

Command entry points exist only for bundled wheels. PYALPS_BUNDLE_APPLICATIONS
is shared with CMake so bindings-only installs cannot overwrite SDK commands.
"""

from __future__ import annotations

import os
import re
import runpy
import sys
from email.parser import BytesParser
from pathlib import Path

# A checkout carries the version file two levels above this project; an
# sdist carries it at the project root through sdist.force-include.
_CORE_PATTERN = r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)"
_CORE = re.compile(_CORE_PATTERN)
_TAG = re.compile(rf"v(?P<core>{_CORE_PATTERN})(?:-(?P<label>(?:alpha|beta|rc|dev)\.[0-9]+))?")

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


def _read_core(project_dir: Path) -> str:
    """Return MAJOR.MINOR.PATCH from ALPS_VERSION.txt."""
    candidates = (project_dir / "ALPS_VERSION.txt", project_dir / "../../cmake/ALPS_VERSION.txt")
    for candidate in candidates:
        if not candidate.is_file():
            continue
        core = candidate.read_text(encoding="utf-8").strip()
        match = _CORE.fullmatch(core)
        if match is None:
            raise RuntimeError(
                f"{candidate} must contain exactly MAJOR.MINOR.PATCH, but reads "
                f"{core!r}. A prerelease label belongs in the "
                f"ALPS_VERSION_PRERELEASE environment variable, and the leading "
                f"'v' of a release tag is not part of the version."
            )
        return core

    tried = ", ".join(str(path) for path in candidates)
    raise RuntimeError(
        "Cannot find ALPS_VERSION.txt, which supplies the pyalps version. "
        f"Looked in: {tried} (relative to {Path.cwd()}). A build from the ALPS "
        "repository finds it under cmake/; an sdist carries a copy at its root, "
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


def version(project_dir: Path | None = None, ref: str | None = None) -> str:
    """The full PEP 440 version for this build."""
    project_dir = project_dir or Path(__file__).resolve().parents[1]
    core = _read_core(project_dir)
    label = os.environ.get("ALPS_VERSION_PRERELEASE", "")
    ref = os.environ.get("GITHUB_REF", "") if ref is None else ref
    if ref.startswith("refs/tags/"):
        tag = ref.removeprefix("refs/tags/")
        match = _TAG.fullmatch(tag)
        if match is None:
            raise ValueError(f"Invalid release tag {tag!r}; expected vMAJOR.MINOR.PATCH[-{{alpha,beta,rc,dev}}.N]")
        if match.group("core") != core:
            raise ValueError(f"Release tag {tag} disagrees with ALPS_VERSION.txt ({core})")
        tag_label = match.group("label") or ""
        if label and _pep440_suffix(label) != _pep440_suffix(tag_label):
            raise ValueError(f"Release tag {tag} disagrees with ALPS_VERSION_PRERELEASE={label!r}")
        label = tag_label
    computed = core + _pep440_suffix(label)
    # Version cannot change between an sdist and wheels rebuilt from it
    # (PEP 643). The build machine's prerelease environment is not present
    # when a user later installs the sdist, so preserve its recorded version.
    pkg_info = project_dir / "PKG-INFO"
    if pkg_info.is_file():
        metadata = BytesParser().parsebytes(pkg_info.read_bytes())
        frozen = metadata.get("Version", "")
        if metadata.get("Name") != "pyalps" or not re.fullmatch(
            re.escape(core) + r"(?:(?:a|b|rc)[0-9]+|\.dev[0-9]+)?", frozen
        ):
            raise ValueError("sdist PKG-INFO disagrees with ALPS_VERSION.txt")
        if (label or ref.startswith("refs/tags/")) and computed != frozen:
            raise ValueError(f"Build version {computed} disagrees with sdist version {frozen}")
        return frozen
    return computed


def bundle_applications() -> bool:
    """One build setting for both CMake's payload and wheel entry points."""
    value = os.environ.get("PYALPS_BUNDLE_APPLICATIONS", "ON").upper()
    if value in {"1", "ON", "YES", "TRUE"}:
        return True
    if value in {"0", "OFF", "NO", "FALSE"}:
        return False
    raise ValueError("PYALPS_BUNDLE_APPLICATIONS must be ON or OFF")


def dynamic_metadata(settings, project):  # noqa: ARG001 - provider protocol
    """scikit-build-core dynamic-metadata 0.3 hook."""
    if settings:
        raise RuntimeError("The alps_version provider takes no settings.")
    scripts = {}
    if bundle_applications():
        cli = runpy.run_path(str(Path(__file__).resolve().parents[1] / "src/pyalps_cli/__init__.py"))
        scripts = {name: "pyalps_cli:main" for name in (*cli["NATIVE_PROGRAMS"], *cli["EXPORTERS"])}
    return {"version": version(), "scripts": scripts}


if __name__ == "__main__":  # a convenience for the release process
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--bundle-applications", action="store_true")
    mode.add_argument("--core", action="store_true", help="print the numeric SDK version")
    args = parser.parse_args()
    print(("ON" if bundle_applications() else "OFF") if args.bundle_applications
          else _read_core(Path(__file__).resolve().parents[1]) if args.core else version())
