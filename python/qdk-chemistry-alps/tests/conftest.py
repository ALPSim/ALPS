# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/conftest.py. Changes: the stack is pyalps
# (models, ed) and qdk-chemistry; the native_dmrg, mps_qdk and qsharp gates
# are replaced by classic_dmrg; the fock_matrix fixture stands in for alpsx's
# Hamiltonian.to_sparse_matrix().
"""Shared fixtures and capability gating.

Every test module here imports this package, which imports pyalps at
module scope, so without the solver stack there is nothing to collect. What
varies between installations that *do* have it is gated by markers,
resolved below:

``pyscf``         needs the pyscf package for a reference value
``classic_dmrg``  needs pyalps built with the classic_dmrg bindings
``local_qdk``     needs a qdk-chemistry providing RdmContainer

There is deliberately no speed tier: every model here is small enough that
the whole suite runs in seconds, and the slowest tests are the plugin
import-order checks, which spawn subprocesses.
"""

from __future__ import annotations

import pytest


def _importable(name: str) -> bool:
    try:
        __import__(name)
    except ImportError:
        return False
    return True


_MISSING = [
    name for name in ("pyalps.models", "pyalps.ed", "qdk_chemistry")
    if not _importable(name)
]

# Drop the solver tests from collection rather than failing to import them.
# (``pytest.importorskip`` does NOT work here: a Skipped raised while pytest
# imports a conftest surfaces as a collection error, not a directory skip.)
# The glob covers ``test_*.py``, so modules added later are handled
# automatically. ``environment_test.py`` is named to fall outside that glob on
# purpose, so the suite still collects something and pytest does not exit with
# "no tests ran".
collect_ignore_glob = ["test_*.py"] if _MISSING else []


def _has(check) -> bool:
    try:
        return bool(check())
    except Exception:
        return False


def _classic_dmrg() -> bool:
    from pyalps import classic_dmrg  # noqa: F401

    return True


def pytest_addoption(parser):
    parser.addoption(
        "--require-solver-stack", action="store_true",
        help="Fail instead of silently skipping validation without pyalps/QDK/classic DMRG",
    )


def pytest_sessionstart(session):
    if session.config.getoption("--require-solver-stack"):
        if _MISSING:
            raise pytest.UsageError(f"missing solver stack: {', '.join(_MISSING)}")
        if not _has(_classic_dmrg):
            raise pytest.UsageError("pyalps must provide the classic_dmrg bindings")


def _rdm_container() -> bool:
    from qdk_chemistry.data import RdmContainer  # noqa: F401

    return True


def _pyscf() -> bool:
    import pyscf  # noqa: F401

    return True


_CAPABILITIES = {
    "classic_dmrg": (
        _classic_dmrg,
        "needs pyalps built with the classic_dmrg bindings",
    ),
    "local_qdk": (
        _rdm_container,
        "needs a qdk-chemistry providing RdmContainer (unreleased as of 2.2.1)",
    ),
    "pyscf": (_pyscf, "needs pyscf"),
}


def pytest_runtest_setup(item: pytest.Item) -> None:
    for name, (check, reason) in _CAPABILITIES.items():
        if name in item.keywords and not _has(check):
            pytest.skip(reason)


# --- shared models ---------------------------------------------------------

@pytest.fixture(scope="session")
def hubbard_qdk():
    """qdk-chemistry Hamiltonian factory for a Hubbard chain."""
    from qdk_chemistry.data import LatticeGraph
    from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

    def make(n_sites: int, u: float = 4.0, t: float = 1.0):
        return create_hubbard_hamiltonian(
            LatticeGraph.chain(n_sites), epsilon=0.0, t=t, U=u
        )

    return make


@pytest.fixture(scope="session")
def hubbard_alps():
    """pyalps Hamiltonian for a Hubbard chain (open boundary)."""
    from pyalps import models

    def make(n_sites: int, u: float = 4.0, t: float = 1.0, mu: float = 0.0):
        return models.hubbard(models.lattice.chain(n_sites, bc="open"), t=t, U=u, mu=mu)

    return make


@pytest.fixture(scope="session")
def fock_matrix():
    """Dense matrix of a pyalps Hamiltonian on its whole product space.

    Rows and columns follow the product order, site 0 the most significant
    local index, whatever order pyalps.ed enumerates the states in. The
    operator need not be Hermitian: it is assembled with operator_matrix on
    an identity problem over the same sites.
    """
    import numpy as np
    from pyalps import ed, models

    def make(H) -> np.ndarray:
        sites = [H.site_for(i) for i in range(H.num_sites)]
        problem = ed.prepare_product(
            models.located(sites[0], "I", 0), sites, [], H.fermion_convention
        )
        states = np.asarray(problem.states)
        index = np.zeros(states.shape[0], dtype=np.int64)
        for site, dim in enumerate(problem.site_dims):
            index = index * int(dim) + states[:, site]
        operator = H.operator_decomposition().op
        matrix = problem.operator_matrix(operator).to_scipy().toarray()
        out = np.zeros_like(matrix)
        out[np.ix_(index, index)] = matrix
        return out

    return make


@pytest.fixture(scope="session")
def integrals():
    """Random 8-fold-symmetric (h1, eri) of a given size."""
    import numpy as np

    def make(n: int, seed: int = 3):
        rng = np.random.default_rng(seed)
        h1 = rng.standard_normal((n, n))
        h1 = 0.5 * (h1 + h1.T)
        e = rng.standard_normal((n, n, n, n)) * 0.2
        # impose the 8-fold symmetry of real chemist integrals
        e = e + e.transpose(1, 0, 2, 3)
        e = e + e.transpose(0, 1, 3, 2)
        e = e + e.transpose(2, 3, 0, 1)
        idx = np.arange(n)
        e[idx, idx, idx, idx] = np.abs(e[idx, idx, idx, idx]) + 1.0
        return h1, e

    return make
