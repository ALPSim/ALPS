# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/src/qdk_chemistry_alps/__init__.py, with the
# AlpsClassicDmrgCalculator export of branch WIP/classic-dmrg-chemistry at
# d4384a7f8. Changes: pyalps instead of alpsx's alps package; the MPS
# AlpsDmrgCalculator is not ported.
"""ALPS solver backends for qdk-chemistry.

Importing this package registers the ALPS algorithms with qdk-chemistry's
registry (the convention released qdk-chemistry uses); on versions with
entry-point discovery, registration happens automatically at
``import qdk_chemistry`` instead and this import is a no-op::

    import qdk_chemistry_alps                       # registers the backends
    from qdk_chemistry.algorithms import create

    calc = create("multi_configuration_calculator", "alps_ed")
    energy, wavefunction = calc.run(hamiltonian, n_alpha, n_beta)

The registered names are ``alps_ed`` (exact diagonalization, returning the
CI vector) and ``alps_classic_dmrg`` (classic block DMRG, returning the CI
vector with the reduced density matrices, or the latter alone for larger
active spaces; registered when pyalps has its ``classic_dmrg`` bindings). :mod:`~qdk_chemistry_alps.interop` converts between qdk-chemistry
and ALPS representations of lattices, spin models and fermionic
Hamiltonians.

Importing this package requires pyalps with its ``models`` and ``ed``
bindings; if they are missing, the import fails with build guidance.
"""

from __future__ import annotations

__version__ = "0.2.0.dev0"

from .plugin import load

# Register first: load() checks for pyalps and raises with build guidance if
# it is missing. The re-exports below pull in modules that import pyalps at
# module scope, so they must not run before that check.
load()

from . import abinitio, interop  # noqa: E402
from .abinitio import build_abinitio_hamiltonian  # noqa: E402
from .classic_dmrg_calculator import AlpsClassicDmrgCalculator  # noqa: E402
from .ed_calculator import AlpsEdCalculator  # noqa: E402

__all__ = [
    "AlpsClassicDmrgCalculator",
    "AlpsEdCalculator",
    "__version__",
    "abinitio",
    "build_abinitio_hamiltonian",
    "interop",
    "load",
]
