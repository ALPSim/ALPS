# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# WIP/classic-dmrg-chemistry at d4384a7f8,
# python/qdk-chemistry-alps/tests/test_plugin.py. Changes: the registered
# names are alps_ed and alps_classic_dmrg (alps_dmrg is not ported); the
# entry-point plugin is also registered through a stand-in registrar.
"""Registration of the ALPS backends with qdk-chemistry.

Two registration paths coexist (import-time and the ``qdk_chemistry.plugins``
entry point) and either may run first depending on which package the user
imports first. Both orders must end with exactly one registration of each
name and no warnings.
"""

from __future__ import annotations

import subprocess
import sys

import pytest

EXPECTED = ["alps_classic_dmrg", "alps_ed"]


def _run(statements: list[str]) -> str:
    """Run statements in a fresh interpreter with warnings as errors.

    A fresh process is the only honest way to test import order: within one
    session both packages are already imported.
    """
    code = "\n".join(
        [
            *statements,
            "from qdk_chemistry.algorithms import available",
            "print(sorted(n for n in available('multi_configuration_calculator')"
            " if n.startswith('alps')))",
        ]
    )
    proc = subprocess.run(
        [sys.executable, "-I", "-W", "error::UserWarning", "-c", code],
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    return proc.stdout


@pytest.mark.classic_dmrg
@pytest.mark.parametrize(
    ("order", "statements"),
    [
        ("alps_first", ["import qdk_chemistry_alps", "import qdk_chemistry"]),
        ("qdk_first", ["import qdk_chemistry", "import qdk_chemistry_alps"]),
    ],
)
def test_registration_is_warning_free_in_either_import_order(order, statements):
    assert _run(statements).strip() == repr(EXPECTED)


@pytest.mark.classic_dmrg
def test_load_is_idempotent():
    from qdk_chemistry.algorithms import available

    import qdk_chemistry_alps

    before = sorted(
        n for n in available("multi_configuration_calculator") if n.startswith("alps")
    )
    for _ in range(3):
        qdk_chemistry_alps.load()
    after = sorted(
        n for n in available("multi_configuration_calculator") if n.startswith("alps")
    )
    assert before == after == EXPECTED


@pytest.mark.classic_dmrg
def test_calculators_resolve_through_the_registry():
    from qdk_chemistry.algorithms import create

    import qdk_chemistry_alps  # noqa: F401  (registers the backends)

    for name in EXPECTED:
        assert create("multi_configuration_calculator", name).name() == name


@pytest.mark.classic_dmrg
def test_classic_calculator_without_rdm_container_says_so(hubbard_qdk):
    """Released qdk-chemistry (2.2.1 and earlier) has no RdmContainer.

    alps_classic_dmrg registers there all the same, so running it has to
    name the missing class rather than fail inside the calculator.
    """
    try:
        from qdk_chemistry.data import RdmContainer  # noqa: F401
    except ImportError:
        pass
    else:
        pytest.skip("this qdk-chemistry provides RdmContainer")
    from qdk_chemistry.algorithms import create

    import qdk_chemistry_alps  # noqa: F401  (registers the backends)

    calculator = create("multi_configuration_calculator", "alps_classic_dmrg")
    with pytest.raises(RuntimeError, match="RdmContainer"):
        calculator.run(hubbard_qdk(2), 1, 1)


def test_entry_point_plugin_class_is_exposed_lazily():
    """AlpsPlugin materializes on attribute access, not at module import.

    The module must not import qdk_chemistry while it is itself being
    imported by qdk-chemistry's plugin discovery.
    """
    pytest.importorskip("qdk_chemistry.plugins", reason="no plugin discovery here")
    from qdk_chemistry_alps import plugin

    assert "AlpsPlugin" in plugin.__all__
    cls = plugin.AlpsPlugin  # triggers the lazy build
    assert cls is plugin.AlpsPlugin  # cached, same object
    assert hasattr(cls, "register")

    with pytest.raises(AttributeError):
        plugin.NoSuchAttribute  # noqa: B018  (attribute access is the test)


def test_entry_point_plugin_skips_registered_names():
    """The entry-point path registers nothing that is already registered."""
    pytest.importorskip("qdk_chemistry.plugins", reason="no plugin discovery here")
    import qdk_chemistry_alps  # noqa: F401  (import-time registration ran)
    from qdk_chemistry_alps import plugin

    class Registrar:
        def __init__(self):
            self.factories = []

        def register_algorithm(self, factory):
            self.factories.append(factory)

    registrar = Registrar()
    plugin.AlpsPlugin().register(registrar)
    assert registrar.factories == []


def test_public_surface():
    import qdk_chemistry_alps as pkg

    for name in pkg.__all__:
        assert hasattr(pkg, name), f"{name} in __all__ but not importable"
    assert pkg.__version__
