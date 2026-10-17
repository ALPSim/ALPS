# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Classic two-site finite-system DMRG on compiled models (``ALPS::classic_dmrg``).

>>> from pyalps import models, classic_dmrg as dmrg
>>> model = models.hubbard(models.lattice.chain(16, bc="open"), t=1.0, U=4.0)
>>> result = dmrg.solve(model, sweeps=dmrg.ramped_sweeps(4, 200), sector={"N": 16, "Sz": 0},
...                     measurements={"n": dmrg.site_terms(model, "n")})
>>> result.energies[0], result.measurements["n"].values[0]      # doctest: +SKIP

``solve`` takes a ``Hamiltonian``, an XML catalog model (``CatalogModel``,
in the sector of its ``*_total`` parameters unless ``sector=`` or ``options``
says otherwise) or a ``CompiledOperator`` with its ``PhysicalBasis``, an
``Options`` object or its fields as keywords, and
measurements as a ``Measurements`` object or a dict of named operators
(measured term by term; use ``expectation=`` for one value per operator).
The operators are anything ``models.compile_observable`` accepts: compiled
operators, located terms or, on a model, operator names (``{"n": "n"}``
measures ``n`` on every site).
Every completed iteration reaches ``callback(iteration)``, which may raise to
abort. The engine is serialized process-wide. The bindings are present when
pyalps was built against an SDK with the ``classic_dmrg`` component.

``keep_mps=True`` also returns every targeted state as a
``MatrixProductState`` in ``result.mps``: the engine's block rotations of
the final sweep and its closing wave function, exactly the state the
energies and measurements belong to. Its ``tensors`` are NumPy arrays
``[left bond, local state, right bond]`` with local states in the
``PhysicalBasis`` order of each site; ``to_dense``, ``amplitude``,
``overlap`` and ``norm`` contract them.

>>> result = dmrg.solve(model, sweeps=dmrg.ramped_sweeps(4, 200), sector={"N": 16, "Sz": 0},
...                     keep_mps=True)                       # doctest: +SKIP
>>> mps = result.mps[0]; mps.bond_dimensions, mps.norm()   # doctest: +SKIP

``evolve(initial, final, time_step=..., steps=...)`` evolves the DMRG ground
state of ``initial`` under ``final`` by adaptive time-dependent DMRG
(White and Feiguin): second-order Trotter gates on the bonds, applied while
sweeping, the kept basis following the state. Real time runs in complex
arithmetic; ``imaginary=True`` applies ``exp(-tau H)`` and renormalizes.
The final Hamiltonian couples only neighbouring sites in chain order;
``drives`` add time-dependent terms. Every measured step reports its time,
the final Hamiltonian's energy, the measurements and, with ``keep_mps``,
its state; ``overlaps`` and ``loschmidt_echo`` contract the kept states.

>>> xx = models.xxz(models.lattice.chain(20, bc="open"), Jxy=1.0, Jz=0.0)
>>> site = xx.site_for(0)
>>> wall = sum(((-1.0 if i < 10 else 1.0) * site.Sz(i) for i in range(20)), models.TermSum())
>>> run = dmrg.evolve(wall, xx, time_step=0.05, steps=60, states=64, sector={"Sz": 0},
...                   measure_steps=10, measurements={"Sz": "Sz"})          # doctest: +SKIP
>>> run.times, run.measurements["Sz"].values                                 # doctest: +SKIP

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

_OPTION_FIELDS = ("sweeps", "warmup_states", "warmup_configuration", "num_eigenvalues", "truncation_error",
                  "truncation_max_states", "lanczos_tolerance", "arithmetic", "sector", "scratch_parent", "verbose",
                  "keep_mps")


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


def make_measurements(measurements=None, expectation=None, model=None):
    """A ``Measurements`` value from dicts of named operators.

    ``measurements`` are reported term by term (one column per term);
    ``expectation`` operators give one value each. Values other than
    compiled operators are compiled on ``model`` (a ``Hamiltonian``,
    ``CatalogModel`` or ``PhysicalBasis``) by ``models.compile_observable``.
    """
    if isinstance(measurements, _classic_dmrg_c.Measurements):
        if expectation:
            raise TypeError("pass either a Measurements object or dicts, not both")
        return measurements

    def compiled(operators):
        out = {}
        for name, value in (operators or {}).items():
            if not isinstance(value, _models.CompiledOperator):
                if model is None:
                    raise TypeError(f"measurement {name!r} needs a model to compile it; pass a CompiledOperator")
                value = _models.compile_observable(model, value)
            out[name] = value
        return out

    result = _classic_dmrg_c.Measurements()
    result.terms = compiled(measurements)
    result.expectation = compiled(expectation)
    return result


def solve(model, basis=None, *, options=None, measurements=None, expectation=None, callback=None,
          log_file=None, **fields):
    """Run the DMRG on a ``Hamiltonian``, a ``CatalogModel`` or a ``CompiledOperator`` with its ``PhysicalBasis``."""
    _models.require_bound(model)
    catalog = getattr(_models, "CatalogModel", None)
    requests = make_measurements(measurements, expectation, model if basis is None else basis)
    if catalog is not None and isinstance(model, catalog):
        if basis is not None:
            raise TypeError("a CatalogModel brings its own physical basis")
        if options is None and "sector" not in fields:
            fields["sector"] = model.sector()
        model, basis = model.compile_hamiltonian(), model.physical_basis()
    resolved = make_options(options, **fields)
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

    ``keep_mps=True`` returns the final state as a ``MatrixProductState``
    over the orbitals in ``result.mps[0]``, on ``spinful_fermion()`` sites
    (local states 0, up, down, up-down) in the same mode order.

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


# -- time evolution -------------------------------------------------------------------

_EVOLUTION_FIELDS = ("time_step", "steps", "imaginary", "states", "truncation_error", "measure_steps", "drives",
                     "keep_mps")


def _is_model(value):
    catalog = getattr(_models, "CatalogModel", None)
    return isinstance(value, _models.Hamiltonian) or (catalog is not None and isinstance(value, catalog))


def _compiled_with_basis(model, basis, role):
    """(CompiledOperator, PhysicalBasis) of a Hamiltonian or a CatalogModel, or of a CompiledOperator or a
    TermSum of located terms on ``basis``."""
    _models.require_bound(model)
    if _is_model(model):
        if isinstance(model, _models.Hamiltonian):
            return _models.compile_hamiltonian(model).op, model.physical_basis()
        return model.compile_hamiltonian(), model.physical_basis()
    if isinstance(model, (_models.CompiledOperator, _models.TermSum)):
        if basis is None:
            raise TypeError(f"the {role} Hamiltonian is an operator: pass a model as the other Hamiltonian "
                            "or its PhysicalBasis as basis=")
        return _models.compile_observable(basis, model), basis
    raise TypeError(f"the {role} Hamiltonian must be a Hamiltonian, a CatalogModel, a CompiledOperator or a "
                    f"TermSum, not {type(model).__name__}")


def _measured_steps(measure_steps, steps):
    import numbers

    if measure_steps is None:
        return []
    if isinstance(measure_steps, numbers.Integral):
        every = int(measure_steps)
        if every < 1:
            raise ValueError("measure_steps as an integer is the measuring interval and must be positive")
        chosen = list(range(0, steps + 1, every))
        if chosen[-1] != steps:
            chosen.append(steps)
        return chosen
    return [int(step) for step in measure_steps]


def _drives(drives, steps, time_step, context):
    out = []
    for item in drives or ():
        if isinstance(item, _classic_dmrg_c.Drive):
            out.append(item)
            continue
        op, coefficients = item
        if _is_model(op):  # a whole model's Hamiltonian, on the evolution's basis
            op, basis = _compiled_with_basis(op, None, "drive")
            if context is not None:
                expected = context if isinstance(context, _models.PhysicalBasis) else context.physical_basis()
                _models.require_same_basis(basis, expected)
        elif not isinstance(op, _models.CompiledOperator):
            op = _models.compile_observable(context, op)
        if callable(coefficients):  # sampled at the step midpoints
            values = [float(coefficients((n + 0.5) * time_step)) for n in range(steps)]
        else:
            values = [float(value) for value in coefficients]
        out.append(_classic_dmrg_c.Drive(op, values))
    return out


def make_evolution(evolution=None, *, context=None, **fields):
    """An ``Evolution`` value: a copy of ``evolution`` (or the defaults) with ``fields`` applied.

    ``measure_steps`` is a list of steps or an interval (every that many
    steps from 0, and the last); ``drives`` are ``Drive`` objects or
    ``(operator, coefficients)`` pairs, the coefficients a sequence with one
    value per step or a function of time sampled at the step midpoints
    ``(n + 1/2) * time_step``. A drive operator is a ``CompiledOperator``,
    anything ``models.compile_observable`` compiles on ``context`` (a model
    or a ``PhysicalBasis``), or a whole model, whose Hamiltonian is added.
    """
    result = _classic_dmrg_c.Evolution()
    if evolution is not None:
        for name in _EVOLUTION_FIELDS:
            setattr(result, name, getattr(evolution, name))
    unknown = set(fields) - set(_EVOLUTION_FIELDS)
    if unknown:
        raise TypeError(f"unknown evolution field {sorted(unknown)[0]!r}; the fields are {', '.join(_EVOLUTION_FIELDS)}")
    for name in ("time_step", "steps", "imaginary", "states", "truncation_error", "keep_mps"):
        if name in fields:
            setattr(result, name, fields[name])
    if "measure_steps" in fields:
        result.measure_steps = _measured_steps(fields["measure_steps"], result.steps)
    if "drives" in fields:
        result.drives = _drives(fields["drives"], result.steps, result.time_step, context)
    return result


def evolve(initial, final, basis=None, *, evolution=None, options=None, measurements=None, expectation=None,
           callback=None, log_file=None, **fields):
    """Time-evolve the DMRG ground state of ``initial`` under ``final``.

    ``initial`` and ``final`` are ``Hamiltonian`` or ``CatalogModel``
    objects on the same physical basis; either may instead be an operator
    (a ``CompiledOperator`` or a ``TermSum`` of located terms) on the
    other's basis or on ``basis``. The initial state is the ground state of ``initial`` in
    the sector (``sector=``, or a catalog model's ``*_total`` parameters),
    solved under ``options`` or the DMRG option fields as for
    :func:`solve`. The evolution fields (``time_step``, ``steps``,
    ``imaginary``, ``states``, ``truncation_error``, ``measure_steps``,
    ``drives``, ``keep_mps``; see :func:`make_evolution`) or an
    ``evolution`` object control the run. ``final`` and the drives must be
    Hermitian, conserve the sector's charges and couple only sites ``i``
    and ``i + 1``. Measurements are named operators as for :func:`solve`.

    Returns an ``EvolutionResult``: ``initial`` (the preparing solve),
    ``steps``, ``times`` and ``energies`` (``<H_final>``) of the measured
    steps, ``measurements`` with one row per measured step,
    ``truncation_errors`` with the discarded weight of every step and, with
    ``keep_mps=True``, the state of every measured step in ``mps``.
    """
    catalog = getattr(_models, "CatalogModel", None)
    both_models = all(isinstance(m, _models.Hamiltonian) for m in (initial, final))
    # Operator names (measurements, drives) resolve on a model where there is one.
    context = next((m for m in (initial, final) if _is_model(m)), None)
    if both_models and basis is None:
        _models.require_bound(initial)
        _models.require_bound(final)
    else:
        if basis is None and context is not None:
            basis = context.physical_basis()
        initial_op, initial_basis = _compiled_with_basis(initial, basis, "initial")
        final_op, final_basis = _compiled_with_basis(final, basis, "final")
        _models.require_same_basis(final_basis, initial_basis)
        if context is None:
            context = initial_basis
    if catalog is not None and isinstance(initial, catalog) and options is None and "sector" not in fields:
        fields["sector"] = initial.sector()
    evolution_fields = {name: fields.pop(name) for name in _EVOLUTION_FIELDS if name in fields}
    if evolution is None and ("time_step" not in evolution_fields or "steps" not in evolution_fields):
        raise TypeError("evolve needs time_step= and steps= (or an Evolution object)")
    resolved_evolution = make_evolution(evolution, context=context, **evolution_fields)
    requests = make_measurements(measurements, expectation, context)
    resolved = make_options(options, **fields)
    if both_models and basis is None:
        return _classic_dmrg_c.evolve(initial, final, resolved_evolution, resolved, requests, callback, log_file)
    return _classic_dmrg_c.evolve(initial_op, final_op, initial_basis, resolved_evolution, resolved, requests,
                                  callback, log_file)


def overlaps(result, reference=None):
    """<reference|psi> for every state an evolution kept (``keep_mps=True``).

    ``reference`` is a ``MatrixProductState`` on the same sites, by default
    the first kept state (the initial state when step 0 is measured).
    """
    import numpy as np

    if not result.mps:
        raise ValueError("the evolution kept no states; run it with keep_mps=True")
    bra = result.mps[0] if reference is None else reference
    return np.array([bra.overlap(state) for state in result.mps], dtype=np.complex128)


def loschmidt_echo(result, reference=None):
    """|<psi(0)|psi(t)>|^2 at every measured step, from the kept states.

    Without ``reference`` step 0 must be measured: its state is psi(0).
    """
    import numpy as np

    if reference is None and (len(result.steps) == 0 or int(result.steps[0]) != 0):
        raise ValueError("the Loschmidt echo needs the initial state: measure step 0 or pass reference=")
    return np.abs(overlaps(result, reference)) ** 2


# -- matrix product states ------------------------------------------------------------

def _mps_to_dense(self, max_dimension=1 << 24):
    """The state on the whole product basis, site 0 most significant (the order of ``pyalps.ed`` states)."""
    import math

    dimension = math.prod(self.local_dimensions)
    if dimension > max_dimension:
        raise ValueError(f"the product basis has {dimension} states, above max_dimension={max_dimension}; "
                         "use amplitude() on chosen product states")
    tensors = self.tensors
    out = tensors[0].reshape(-1, tensors[0].shape[2])
    for tensor in tensors[1:]:
        out = (out @ tensor.reshape(tensor.shape[0], -1)).reshape(-1, tensor.shape[2])
    return out.reshape(-1)


def _mps_amplitude(self, states):
    """<s|psi> for one product state (a sequence of local indices) or for each row of an (n, L) array."""
    import numpy as np

    configurations = np.asarray(states, dtype=np.int64)
    single = configurations.ndim == 1
    configurations = np.atleast_2d(configurations)
    tensors = self.tensors
    if configurations.shape[1] != len(tensors):
        raise ValueError(f"a product state of this MPS has {len(tensors)} local indices")
    carry = np.ones((configurations.shape[0], 1), dtype=np.result_type(*tensors))
    for site, tensor in enumerate(tensors):
        local = configurations[:, site]
        if local.min(initial=0) < 0 or local.max(initial=0) >= tensor.shape[1]:
            raise ValueError(f"a local index at site {site} is outside 0..{tensor.shape[1] - 1}")
        carry = np.einsum("nl,lnr->nr", carry, tensor[:, local, :])
    values = carry[:, 0]
    return values[0] if single else values


def _mps_overlap(self, other):
    """<self|other>, contracted site by site; both must have the same local dimensions."""
    import numpy as np

    if list(self.local_dimensions) != list(other.local_dimensions):
        raise ValueError("the matrix product states have different local dimensions")
    environment = np.ones((1, 1))
    for a, b in zip(self.tensors, other.tensors):
        environment = np.einsum("xy,xsl,ysr->lr", environment, a.conj(), b)
    return environment[0, 0]


def _mps_norm(self):
    """sqrt(<psi|psi>); one for a state from the solver."""
    import numpy as np

    return float(np.sqrt(abs(_mps_overlap(self, self))))


_classic_dmrg_c.MatrixProductState.to_dense = _mps_to_dense
_classic_dmrg_c.MatrixProductState.amplitude = _mps_amplitude
_classic_dmrg_c.MatrixProductState.overlap = _mps_overlap
_classic_dmrg_c.MatrixProductState.norm = _mps_norm


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
    "make_options", "make_measurements", "solve", "make_evolution", "evolve", "overlaps", "loschmidt_echo",
    "local_expectations", "correlation_matrix", "spin_spin_correlations", "ci_occupations",
]
