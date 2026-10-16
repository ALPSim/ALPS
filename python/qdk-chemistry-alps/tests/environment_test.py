# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/environment_test.py. Changes: the stack is
# pyalps; the documentation build check is gone with the mkdocs site, which
# is not ported.
"""Checks that run with or without the solver stack installed.

Named ``environment_test.py`` rather than ``test_environment.py`` on
purpose: ``conftest.py`` drops ``test_*.py`` from collection when pyalps
or ``qdk_chemistry`` is missing, and this module has to survive that. It
keeps the suite non-empty in a bare environment (so pytest does not exit
"no tests ran") and it verifies the thing that matters most to someone
installing for the first time — that a missing dependency says so clearly.
"""

from __future__ import annotations

import importlib.metadata
import subprocess
import sys

import pytest


def _importable(name: str) -> bool:
    # Duplicated from conftest rather than imported: tests/ is not a package,
    # so a relative import here would depend on rootdir layout.
    try:
        __import__(name)
    except ImportError:
        return False
    return True


HAS_STACK = (
    _importable("pyalps.models") and _importable("pyalps.ed")
    and _importable("qdk_chemistry")
)


def test_distribution_metadata_is_installed():
    """Readable without importing the package, so it works in any environment."""
    metadata = importlib.metadata.metadata("qdk-chemistry-alps")
    assert metadata["Name"] == "qdk-chemistry-alps"
    assert importlib.metadata.version("qdk-chemistry-alps")


def test_entry_point_is_declared():
    points = importlib.metadata.entry_points(group="qdk_chemistry.plugins")
    assert any(point.name == "alps" for point in points)


def test_missing_pyalps_is_reported_with_guidance():
    """The first thing a new user hits must point at the fix, not a traceback.

    pyalps is hidden in a fresh interpreter, so this runs where it is
    installed too.
    """
    code = (
        "import sys\n"
        "sys.modules['pyalps'] = None\n"
        "try:\n"
        "    import qdk_chemistry_alps\n"
        "except ImportError as exc:\n"
        "    print(exc)\n"
        "else:\n"
        "    raise SystemExit('imported without pyalps')\n"
    )
    proc = subprocess.run(
        [sys.executable, "-I", "-c", code], capture_output=True, text=True, timeout=300
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "pyalps" in proc.stdout
    assert "README" in proc.stdout


@pytest.mark.skipif(not HAS_STACK, reason="needs pyalps and qdk_chemistry")
def test_backends_register_when_the_stack_is_present():
    from qdk_chemistry.algorithms import available

    import qdk_chemistry_alps  # noqa: F401

    names = set(available("multi_configuration_calculator"))
    assert "alps_ed" in names
    assert ("alps_classic_dmrg" in names) == _importable("pyalps.classic_dmrg")
