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
matrices in the PySCF convention, their spin blocks, any further
``observables`` and, on request, the CI vector. ``local_expectations``,
``correlation_matrix`` and ``spin_spin_correlations`` derive orbital-resolved
quantities from those RDMs; ``ci_occupations`` unpacks the CI determinants.
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


def solve_molecular(one_body, two_body, n_alpha, n_beta, *, core_energy=0.0, integral_cutoff=0.0, observables=None,
                    ci_vector=False, options=None, callback=None, log_file=None, **fields):
    """Classic block DMRG of a restricted real molecular Hamiltonian.

    ``one_body[p, q]`` and ``two_body[p, q, r, s]`` are spatial-orbital
    integrals in chemist notation, ``(pq|rs)``, with the real-orbital
    eightfold symmetry after ``integral_cutoff`` drops smaller entries. The
    ground state of the ``(n_alpha, n_beta)`` sector is returned as a
    ``MolecularResult``: ``energy`` (including ``core_energy``), the spin
    blocks ``one_rdm_alpha``/``one_rdm_beta`` with ``dm1[p, q] = <c+_q c_p>``
    and ``two_rdm[p, q, r, s] = sum_ab <c+_pa c+_rb c_sb c_qa>``, with its
    spin blocks ``two_rdm_aa``/``two_rdm_ab``/``two_rdm_bb`` as PySCF's
    ``make_rdm12s``. ``options`` or keyword fields set the sweeps, warm-up,
    truncation and Lanczos controls as for :func:`solve`; ``sector``,
    ``num_eigenvalues`` and ``arithmetic`` must keep their defaults. RDM
    storage grows as the fourth power of the orbital count.

    ``observables`` maps names to ``pyalps.models`` term sums on
    ``spinful_fermion()`` sites, site ``p`` being orbital ``p``, with real
    coefficients and at most four sites per product; they are measured on
    the final state with the RDMs and returned in ``result.observables``.
    Nothing can be measured after the call, as the engine keeps no state.

    ``ci_vector=True`` also returns the final state on every determinant of
    the sector (``ci_alpha``, ``ci_beta``, ``ci_coefficients``; see
    :func:`ci_occupations`), contracted from the engine's stored rotations.
    The amplitudes use the ALPS mode order of ``pyalps.ed``: orbitals
    ascending, each orbital's down mode before its up mode. The sector grows
    combinatorially, so this suits small active spaces.
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
                                           dict(observables or {}), bool(ci_vector), resolved, callback, log_file)


# Orbital-resolved quantities of a MolecularResult, derived from its RDMs.
# The number-type operators are combinations (up, down) of n_up and n_dn.
_NUMBER_OPERATORS = {"n": (1.0, 1.0), "n_up": (1.0, 0.0), "n_dn": (0.0, 1.0), "Sz": (0.5, -0.5)}
_LADDER_OPERATORS = ("c_up", "c_dn", "Sp", "Sm")


def _number_pairs(result):
    """<n[i, s] n[j, t]> for s, t in (up, down), from the spin-resolved RDMs."""
    import numpy as np

    a, b = np.diag(result.one_rdm_alpha), np.diag(result.one_rdm_beta)
    aa, ab, bb = result.two_rdm_aa, result.two_rdm_ab, result.two_rdm_bb
    up_up = np.einsum("iijj->ij", aa) + np.diag(a)
    down_down = np.einsum("iijj->ij", bb) + np.diag(b)
    up_down = np.einsum("iijj->ij", ab)
    return {("up", "up"): up_up, ("up", "down"): up_down, ("down", "up"): up_down.T, ("down", "down"): down_down}


def local_expectations(result, op):
    """<O_i> for every orbital: ``op`` is one of n, n_up, n_dn, Sz or n_updn."""
    import numpy as np

    if op == "n_updn":
        return np.einsum("iiii->i", result.two_rdm_ab).copy()
    if op not in _NUMBER_OPERATORS:
        raise ValueError(f"local_expectations knows n, n_up, n_dn, Sz and n_updn, not {op!r}; "
                         "measure other operators through solve_molecular(observables=...)")
    up, down = _NUMBER_OPERATORS[op]
    return up * np.diag(result.one_rdm_alpha) + down * np.diag(result.one_rdm_beta)


def correlation_matrix(result, op_i, op_j=None):
    """C[i, j] = <O_i+ P_j> over all orbital pairs, from the RDMs.

    ``op_i`` and ``op_j`` (default ``op_i``) are number-type operators (n,
    n_up, n_dn, Sz), in any combination, or one ladder operator (c_up, c_dn,
    Sp, Sm) on both sides, its Jordan-Wigner string implied:
    ``correlation_matrix(r, "c_up")`` is ``<c+[i, up] c[j, up]>`` and
    ``correlation_matrix(r, "Sp")`` is ``<Sm_i Sp_j>``.
    """
    import numpy as np

    op_j = op_i if op_j is None else op_j
    if op_i in _NUMBER_OPERATORS and op_j in _NUMBER_OPERATORS:
        pairs = _number_pairs(result)
        left, right = dict(zip(("up", "down"), _NUMBER_OPERATORS[op_i])), dict(zip(("up", "down"), _NUMBER_OPERATORS[op_j]))
        return sum(left[s] * right[t] * pairs[(s, t)] for s in ("up", "down") for t in ("up", "down"))
    if op_i == op_j and op_i in _LADDER_OPERATORS:
        alpha, beta, ab = result.one_rdm_alpha, result.one_rdm_beta, result.two_rdm_ab
        if op_i == "c_up":
            return alpha.T.copy()
        if op_i == "c_dn":
            return beta.T.copy()
        if op_i == "Sp":  # <Sm_i Sp_j> = d_ij <n_dn> - ab[j, i, i, j]
            return np.diag(np.diag(beta)) - np.einsum("jiij->ij", ab)
        return np.diag(np.diag(alpha)) - np.einsum("ijji->ij", ab)  # Sm: <Sp_i Sm_j>
    raise ValueError(f"correlation_matrix cannot derive <{op_i}+ {op_j}> from the RDMs; "
                     "measure it through solve_molecular(observables=...)")


def spin_spin_correlations(result):
    """<S_i . S_j> over all orbital pairs; its sum is <S^2>."""
    return correlation_matrix(result, "Sz") + 0.5 * (correlation_matrix(result, "Sp") + correlation_matrix(result, "Sm"))


def ci_occupations(result):
    """The CI determinants as two (determinants, orbitals) 0/1 arrays, alpha and beta."""
    import numpy as np

    bits = np.arange(result.orbitals, dtype=np.uint64)
    alpha = ((result.ci_alpha[:, None] >> bits) & np.uint64(1)).astype(np.int8)
    beta = ((result.ci_beta[:, None] >> bits) & np.uint64(1)).astype(np.int8)
    return alpha, beta


__all__ = [name for name in dir(_classic_dmrg_c) if not name.startswith("_")] + [
    "make_options", "make_measurements", "solve", "local_expectations", "correlation_matrix",
    "spin_spin_correlations", "ci_occupations",
]
