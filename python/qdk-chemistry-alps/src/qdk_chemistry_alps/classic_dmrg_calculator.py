# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# WIP/classic-dmrg-chemistry at d4384a7f8,
# python/qdk-chemistry-alps/src/qdk_chemistry_alps/classic_dmrg_calculator.py.
# Changes: the solve calls pyalps.classic_dmrg.solve_molecular (alpsx:
# alps.solvers.classic_dmrg.run_molecular, the same C++ solver); the
# integrals are passed as contiguous float64 arrays.
#
# The wavefunction_container, ci_coefficient_cutoff and ci_memory_limit_gb
# settings, the container dispatch, the memory guard, the determinant loop
# and measure_local / measure_correlations / spin_spin_correlations are
# adapted from alpsx's MPS calculator, branch agent/v3-execution-and-api-audit
# at d37e51227, python/qdk-chemistry-alps/src/qdk_chemistry_alps/
# dmrg_calculator.py. Changes: the container is chosen before the solve,
# because the CI vector is requested from solve_molecular, which contracts it
# from the engine's stored rotations (alpsx contracted a retained MPS
# afterwards); there is no "mps" mode; the memory estimate is per
# determinant; the measurements are derived from the run's RDMs (alpsx
# contracted its MPS), and other operators are measured in the final sweep
# through the observables attribute; the run's spin-resolved RDMs are
# attached to both containers.
"""QDK molecular calculator using the classic ALPS block DMRG engine."""

from __future__ import annotations

from math import comb

import numpy as np
from qdk_chemistry.algorithms import MultiConfigurationCalculator
from qdk_chemistry.data import (
    Configuration,
    Settings,
    StateVectorContainer,
    Wavefunction,
)

from .abinitio import interleaved_to_blocked_parity

__all__ = ["AlpsClassicDmrgCalculator"]

# Spin-half string characters by alpha + 2 * beta occupation.
_STATE_CHAR = np.array(["0", "u", "d", "2"])
_MODES = ("auto", "statevector", "rdm")
# Peak memory of a CI vector, per determinant: about 460 bytes measured for
# the QDK Configuration objects and container, plus the engine's
# reconstruction and result arrays and their NumPy copies, rounded up.
_CI_BYTES_PER_DETERMINANT = 1000
# The engine packs a determinant's occupations into 64-bit words.
_CI_MAX_ORBITALS = 64


class AlpsClassicDmrgCalculator(MultiConfigurationCalculator):
    """Restricted real molecular DMRG through pyalps's classic block engine.

    pyalps's ``classic_dmrg.solve_molecular`` builds the molecular
    Hamiltonian, runs the sweeps and measures, in the final sweep, the 1-/2-RDMs
    by spin block and, on request, the CI vector of the final state,
    contracted from the engine's stored rotations. The engine keeps no
    state afterwards. The container depends on ``wavefunction_container``:

    * ``"auto"`` (default): a ``StateVectorContainer`` while the sector's CI
      vector fits in ``ci_memory_limit_gb``, otherwise an ``RdmContainer``.
    * ``"statevector"``: always the CI vector; ``MemoryError`` beyond the
      limit.
    * ``"rdm"``: an ``RdmContainer`` without a CI vector.

    ``RdmContainer`` needs a qdk-chemistry build that provides it (released
    2.2.1 does not); the state-vector path works with released builds. The
    ``StateVectorContainer`` carries the run's RDMs unless
    ``ci_coefficient_cutoff`` trims the vector, whose RDMs then differ.
    The engine keeps complementary block operators, but the CI vector grows
    combinatorially: at 12 orbitals and half filling, building the
    ``StateVectorContainer`` takes about as long as the solve, and ``"rdm"``
    avoids it.

    :meth:`make_rdms`, :meth:`measure_local`, :meth:`measure_correlations`
    and :meth:`spin_spin_correlations` read the last run. Other operators are
    measured in the final sweep: set :attr:`observables` to a dict of
    ``pyalps.models`` term sums on ``spinful_fermion()`` sites (site ``p`` is
    orbital ``p``) before :meth:`run`, and read the values from
    ``last_result.observables``.
    """

    def __init__(self):
        super().__init__()
        self._settings = Settings()
        self._settings._set_default("integral_cutoff", "double", 1e-12)
        self._settings._set_default("max_bond_dim", "int", 200)
        self._settings._set_default("warmup_states", "int", 20)
        self._settings._set_default("sweeps", "int", 20)
        self._settings._set_default("convergence_tol", "double", 1e-10)
        self._settings._set_default("truncation_tol", "double", 0.0)
        self._settings._set_default("ci_coefficient_cutoff", "double", 0.0)
        self._settings._set_default("ci_memory_limit_gb", "double", 4.0)
        self._settings._set_default("wavefunction_container", "string", "auto")
        self.observables = {}
        self.last_result = None
        self.last_energies = None

    def name(self) -> str:
        return "alps_classic_dmrg"

    def _last(self):
        if self.last_result is None:
            raise RuntimeError("run() has not produced a classic DMRG result yet")
        return self.last_result

    def make_rdms(self) -> tuple[np.ndarray, np.ndarray]:
        """Detached spin-summed 1-/2-RDMs in the PySCF convention."""
        result = self._last()
        return result.one_rdm, result.two_rdm.copy()

    def measure_local(self, op: str) -> np.ndarray:
        """Orbital-resolved <O_i> of the last run, from its RDMs.

        ``op`` is ``"n"``, ``"n_up"``, ``"n_dn"``, ``"Sz"`` or ``"n_updn"``
        (double occupancy).
        """
        from pyalps import classic_dmrg

        return self._derived(classic_dmrg.local_expectations, op)

    def measure_correlations(self, op_i: str, op_j: str | None = None) -> np.ndarray:
        """Two-point correlation matrix C_ij = <O_i+ O_j> of the last run.

        Derived from its RDMs: ``op_i`` and ``op_j`` (default ``op_i``) are
        number-type operators (``"n"``, ``"n_up"``, ``"n_dn"``, ``"Sz"``) in
        any combination, or one of ``"c_up"``, ``"c_dn"``, ``"Sp"``, ``"Sm"``
        on both sides, Jordan-Wigner strings implied::

            calc.measure_correlations("Sz")    # <Sz_i Sz_j>
            calc.measure_correlations("n")     # <n_i n_j>
            calc.measure_correlations("c_up")  # <c+_{i,up} c_{j,up}>
            calc.measure_correlations("Sp")    # <Sm_i Sp_j>, as (Sp)+ = Sm
        """
        from pyalps import classic_dmrg

        return self._derived(classic_dmrg.correlation_matrix, op_i, op_j)

    def spin_spin_correlations(self) -> np.ndarray:
        """<S_i . S_j> over all orbital pairs of the last run; its sum is <S^2>."""
        from pyalps import classic_dmrg

        return classic_dmrg.spin_spin_correlations(self._last())

    def _derived(self, function, *ops):
        result = self._last()
        try:
            return function(result, *ops)
        except ValueError as exc:
            raise ValueError(
                f"{exc}; with alps_classic_dmrg, set calculator.observables before "
                "run() and read last_result.observables"
            ) from None

    def _run_impl(self, hamiltonian, n_alpha: int, n_beta: int):
        from pyalps import classic_dmrg

        if not hamiltonian.is_restricted():
            raise NotImplementedError("alps_classic_dmrg supports restricted Hamiltonians only")
        one = np.asarray(hamiltonian.get_one_body_integrals()[0], dtype=np.float64)
        two = hamiltonian.get_two_body_integrals()
        two = two[0] if isinstance(two, (tuple, list)) else two
        n = one.shape[0]
        statevector = self._wants_statevector(n, n_alpha, n_beta)
        if not statevector:
            _rdm_container_class()  # fail before the solve, not after it
        truncation = float(self._settings.get("truncation_tol"))
        # Zero selects the classic fixed-state schedule. Positive values opt
        # into error targeting, which can reduce the retained block spaces.
        error_target = None if truncation == 0.0 else truncation
        max_states = int(self._settings.get("max_bond_dim"))
        result = classic_dmrg.solve_molecular(
            np.ascontiguousarray(one),
            np.ascontiguousarray(np.asarray(two, dtype=np.float64).reshape(n, n, n, n)),
            n_alpha,
            n_beta,
            core_energy=float(hamiltonian.get_core_energy()),
            integral_cutoff=float(self._settings.get("integral_cutoff")),
            observables=dict(self.observables),
            ci_vector=statevector,
            sweeps=classic_dmrg.fixed_sweeps(int(self._settings.get("sweeps")), max_states),
            warmup_states=int(self._settings.get("warmup_states")),
            lanczos_tolerance=float(self._settings.get("convergence_tol")),
            truncation_error=error_target,
            truncation_max_states=None if error_target is None else max_states,
        )
        self.last_result = result
        self.last_energies = np.asarray([result.energy])
        if statevector:
            container = self._statevector_container(hamiltonian, result)
        else:
            container = self._rdm_container(hamiltonian, result, n_alpha, n_beta, max_states)
        return result.energy, Wavefunction(container)

    def _ci_workspace_gb(self, n, n_alpha, n_beta) -> float:
        return comb(n, n_alpha) * comb(n, n_beta) * _CI_BYTES_PER_DETERMINANT / 1e9

    def _wants_statevector(self, n, n_alpha, n_beta) -> bool:
        """Dispatch on ``wavefunction_container``, before the solve."""
        mode = str(self._settings.get("wavefunction_container"))
        if mode not in _MODES:
            raise ValueError(
                f"wavefunction_container={mode!r}; expected 'auto', 'statevector' or 'rdm'"
            )
        if mode == "rdm":
            return False
        workspace_gb = self._ci_workspace_gb(n, n_alpha, n_beta)
        limit = float(self._settings.get("ci_memory_limit_gb"))
        fits = n <= _CI_MAX_ORBITALS and workspace_gb <= limit
        if mode == "auto" or fits:
            return fits
        if n > _CI_MAX_ORBITALS:
            raise ValueError(
                f"the classic CI reconstruction handles at most {_CI_MAX_ORBITALS} "
                f"orbitals, not {n}; use wavefunction_container='rdm' or 'auto'"
            )
        raise MemoryError(
            f"the CI vector of the ({n_alpha},{n_beta}) sector on {n} orbitals has "
            f"{comb(n, n_alpha) * comb(n, n_beta)} determinants (~{workspace_gb:.1f} GB > "
            f"ci_memory_limit_gb={limit}); use wavefunction_container='rdm' (or "
            "'auto') for large active spaces, or raise the limit"
        )

    def _rdm_data(self, result):
        """The run's RDMs as QDK container arguments.

        QDK's aaaa/aabb/bbbb blocks are pyalps's two_rdm_aa/ab/bb as they
        are: Gamma_aabb[p, q, r, s] = <c+_{p,a} c+_{r,b} c_{s,b} c_{q,a}>
        (checked element by element against MACIS's spin-dependent 2-RDM).
        The spin-resolved 1-RDMs matter because model orbitals carry no spin
        axis, so a container cannot split the spin-traced one.
        """
        def flat(array):
            return np.ascontiguousarray(array).reshape(-1)

        return dict(
            one_rdm_spin_traced=np.ascontiguousarray(result.one_rdm),
            one_rdm_aa=np.ascontiguousarray(result.one_rdm_alpha),
            one_rdm_bb=np.ascontiguousarray(result.one_rdm_beta),
            two_rdm_spin_traced=flat(result.two_rdm),
            two_rdm_aaaa=flat(result.two_rdm_aa),
            two_rdm_aabb=flat(result.two_rdm_ab),
            two_rdm_bbbb=flat(result.two_rdm_bb),
        )

    def _rdm_container(self, hamiltonian, result, n_alpha, n_beta, max_states):
        metadata = {
            "energy": result.energy,
            "max_bond_dim": max_states,
            "sweeps": int(self._settings.get("sweeps")),
            "final_truncation_error": result.final_truncation_error,
        }
        return _rdm_container_class()(
            hamiltonian.get_orbitals(), n_alpha, n_beta,
            **self._rdm_data(result), method=self.name(), metadata=metadata,
        )

    def _statevector_container(self, hamiltonian, result):
        """The run's CI vector, in QDK's determinant order and signs."""
        from pyalps import classic_dmrg

        # Changed from alpsx's per-determinant loop: the signs and the
        # spin-half strings are formed for all kept determinants at once,
        # which took most of the time at 12 orbitals (853776 determinants).
        alpha, beta = classic_dmrg.ci_occupations(result)
        amplitudes = np.asarray(result.ci_coefficients)
        cutoff = float(self._settings.get("ci_coefficient_cutoff"))
        order = np.argsort(-np.abs(amplitudes), kind="stable")
        if cutoff > 0.0:
            order = order[np.abs(amplitudes[order]) >= cutoff]
        alpha, beta = alpha[order], beta[order]
        coefficients = amplitudes[order] * interleaved_to_blocked_parity(alpha, beta)
        # Joining a row of one-character strings is a view as one string per row.
        strings = np.ascontiguousarray(_STATE_CHAR[alpha + 2 * beta]).view(f"<U{alpha.shape[1]}").ravel()
        determinants = [Configuration.from_spin_half_string(s) for s in strings.tolist()]
        trimmed = len(coefficients) < len(amplitudes)
        coefficients /= np.linalg.norm(coefficients)
        if trimmed:
            return StateVectorContainer(coefficients, determinants, hamiltonian.get_orbitals())
        return StateVectorContainer(
            coefficients, determinants, hamiltonian.get_orbitals(), **self._rdm_data(result)
        )


def _rdm_container_class():
    try:
        from qdk_chemistry.data import RdmContainer
    except ImportError as exc:
        raise RuntimeError(
            "this qdk-chemistry build has no RdmContainer, which alps_classic_dmrg "
            "needs without a CI vector; use wavefunction_container='statevector' "
            "(or 'auto' with a larger ci_memory_limit_gb) for small active spaces, "
            "or pyalps.classic_dmrg.solve_molecular for the RDMs alone"
        ) from exc
    return RdmContainer
