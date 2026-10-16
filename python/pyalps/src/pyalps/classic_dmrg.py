# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Classic two-site finite-system DMRG on compiled models (``ALPS::classic_dmrg``).

>>> from pyalps import models, classic_dmrg as dmrg
>>> model = models.hubbard(models.lattice.chain(16, bc="open"), t=1.0, U=4.0)
>>> result = dmrg.solve(model, sweeps=dmrg.ramped_sweeps(4, 200), sector={"N": 16, "Sz": 0},
...                     measurements={"n": dmrg.site_terms(model, "n")})
>>> result.energies[0], result.measurements["n"].values[0]      # doctest: +SKIP

``solve`` takes a ``Hamiltonian`` (or a ``CompiledOperator`` with its
``PhysicalBasis``), an ``Options`` object or its fields as keywords, and
measurements as a ``Measurements`` object or a dict of named operators
(measured term by term; use ``expectation=`` for one value per operator).
Every completed iteration reaches ``callback(iteration)``, which may raise to
abort. The engine is serialized process-wide. The bindings are present when
pyalps was built against an SDK with the ``classic_dmrg`` component.

``solve_molecular(h1, eri, n_alpha, n_beta)`` runs the same engine on a
restricted real molecular Hamiltonian (chemist-notation integrals) and
returns the energy with the final-state one- and two-particle density
matrices in the PySCF convention.
"""

from . import models as _models  # registers the model types first  # noqa: F401

try:
    from ._ext import classic_dmrg_c as _classic_dmrg_c
except ImportError as error:  # pragma: no cover - SDK without classic_dmrg
    raise ImportError(
        "pyalps was built without the classic DMRG bindings (PYALPS_BUILD_MODELS=OFF "
        "or an SDK without ALPS::classic_dmrg)"
    ) from error

from ._ext.classic_dmrg_c import *  # noqa: F401,F403

_OPTION_FIELDS = ("sweeps", "warmup_states", "num_eigenvalues", "truncation_error", "truncation_max_states",
                  "lanczos_tolerance", "arithmetic", "sector", "scratch_parent", "verbose")


def make_options(options=None, **fields):
    """An ``Options`` value: a copy of ``options`` (or the defaults) with ``fields`` applied."""
    result = _classic_dmrg_c.Options()
    if options is not None:
        for name in _OPTION_FIELDS:
            setattr(result, name, getattr(options, name))
    for name, value in fields.items():
        if name not in _OPTION_FIELDS:
            raise TypeError(f"unknown DMRG option {name!r}; the options are {', '.join(_OPTION_FIELDS)}")
        setattr(result, name, value)
    return result


def make_measurements(measurements=None, expectation=None):
    """A ``Measurements`` value from dicts of named compiled operators.

    ``measurements`` are reported term by term (one column per term);
    ``expectation`` operators give one value each.
    """
    if isinstance(measurements, _classic_dmrg_c.Measurements):
        if expectation:
            raise TypeError("pass either a Measurements object or dicts, not both")
        return measurements
    result = _classic_dmrg_c.Measurements()
    result.terms = dict(measurements or {})
    result.expectation = dict(expectation or {})
    return result


def solve(model, basis=None, *, options=None, measurements=None, expectation=None, callback=None,
          log_file=None, **fields):
    """Run the DMRG on a ``Hamiltonian`` or on a ``CompiledOperator`` with its ``PhysicalBasis``."""
    resolved = make_options(options, **fields)
    requests = make_measurements(measurements, expectation)
    if basis is None:
        return _classic_dmrg_c.solve(model, resolved, requests, callback, log_file)
    return _classic_dmrg_c.solve(model, basis, resolved, requests, callback, log_file)


def solve_molecular(one_body, two_body, n_alpha, n_beta, *, core_energy=0.0, integral_cutoff=0.0, options=None,
                    callback=None, log_file=None, **fields):
    """Classic block DMRG of a restricted real molecular Hamiltonian.

    ``one_body[p, q]`` and ``two_body[p, q, r, s]`` are spatial-orbital
    integrals in chemist notation, ``(pq|rs)``, with the real-orbital
    eightfold symmetry after ``integral_cutoff`` drops smaller entries. The
    ground state of the ``(n_alpha, n_beta)`` sector is returned as a
    ``MolecularResult``: ``energy`` (including ``core_energy``), the spin
    blocks ``one_rdm_alpha``/``one_rdm_beta`` with ``dm1[p, q] = <c+_q c_p>``
    and ``two_rdm[p, q, r, s] = sum_ab <c+_pa c+_rb c_sb c_qa>``. ``options``
    or keyword fields set the sweeps, warm-up, truncation and Lanczos
    controls as for :func:`solve`; ``sector``, ``num_eigenvalues`` and
    ``arithmetic`` must keep their defaults. RDM storage grows as the fourth
    power of the orbital count.
    """
    import numpy as np

    arrays = []
    for value, name, rank in ((one_body, "one_body", 2), (two_body, "two_body", 4)):
        value = np.asarray(value)
        if value.ndim != rank:
            raise ValueError(f"{name} must have {rank} dimensions")
        if value.dtype.kind not in "fiu":
            raise TypeError(f"{name} must hold real integrals")
        arrays.append(np.ascontiguousarray(value, dtype=np.float64))
    resolved = make_options(options, **fields)
    return _classic_dmrg_c.solve_molecular(*arrays, n_alpha, n_beta, float(core_energy), float(integral_cutoff),
                                           resolved, callback, log_file)


__all__ = [name for name in dir(_classic_dmrg_c) if not name.startswith("_")] + [
    "make_options", "make_measurements", "solve",
]
