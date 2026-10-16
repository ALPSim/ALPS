# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# WIP/classic-dmrg-chemistry at d4384a7f8,
# python/qdk-chemistry-alps/src/qdk_chemistry_alps/classic_dmrg_calculator.py.
# Changes: the solve calls pyalps.classic_dmrg.solve_molecular (alpsx:
# alps.solvers.classic_dmrg.run_molecular, the same C++ solver); the
# integrals are passed as contiguous float64 arrays.
"""QDK molecular calculator using the classic ALPS block DMRG engine."""

from __future__ import annotations

import numpy as np
from qdk_chemistry.algorithms import MultiConfigurationCalculator
from qdk_chemistry.data import Settings, Wavefunction

__all__ = ["AlpsClassicDmrgCalculator"]


class AlpsClassicDmrgCalculator(MultiConfigurationCalculator):
    """Restricted real molecular DMRG returning an RDM-backed wavefunction.

    Requires a qdk-chemistry build providing ``RdmContainer``. pyalps's
    ``classic_dmrg.solve_molecular`` builds the molecular Hamiltonian, runs
    the sweeps and measures the reduced density matrices in C++. The result
    retains 1-/2-RDMs and diagnostics; it does not retain an MPS or
    reconstruct determinant coefficients. Every integral is a separate
    engine term, so the cost grows quickly with the orbital count: use it for
    small active spaces.
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
        self.last_result = None
        self.last_energies = None

    def name(self) -> str:
        return "alps_classic_dmrg"

    def make_rdms(self) -> tuple[np.ndarray, np.ndarray]:
        """Detached spin-summed 1-/2-RDMs in the PySCF convention."""
        if self.last_result is None:
            raise RuntimeError("run() has not produced a classic DMRG result yet")
        return self.last_result.one_rdm, self.last_result.two_rdm.copy()

    def _run_impl(self, hamiltonian, n_alpha: int, n_beta: int):
        from pyalps import classic_dmrg

        try:
            from qdk_chemistry.data import RdmContainer
        except ImportError as exc:
            raise RuntimeError(
                "alps_classic_dmrg requires a qdk-chemistry build with RdmContainer"
            ) from exc
        if not hamiltonian.is_restricted():
            raise NotImplementedError("alps_classic_dmrg supports restricted Hamiltonians only")
        one = np.asarray(hamiltonian.get_one_body_integrals()[0], dtype=np.float64)
        two = hamiltonian.get_two_body_integrals()
        two = two[0] if isinstance(two, (tuple, list)) else two
        n = one.shape[0]
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
            sweeps=classic_dmrg.fixed_sweeps(int(self._settings.get("sweeps")), max_states),
            warmup_states=int(self._settings.get("warmup_states")),
            lanczos_tolerance=float(self._settings.get("convergence_tol")),
            truncation_error=error_target,
            truncation_max_states=None if error_target is None else max_states,
        )
        metadata = {
            "energy": result.energy,
            "max_bond_dim": max_states,
            "sweeps": int(self._settings.get("sweeps")),
            "final_truncation_error": result.final_truncation_error,
        }
        container = RdmContainer(
            hamiltonian.get_orbitals(), n_alpha, n_beta,
            one_rdm_spin_traced=np.ascontiguousarray(result.one_rdm),
            one_rdm_aa=np.ascontiguousarray(result.one_rdm_alpha),
            one_rdm_bb=np.ascontiguousarray(result.one_rdm_beta),
            two_rdm_spin_traced=np.ascontiguousarray(result.two_rdm).reshape(-1),
            method=self.name(), metadata=metadata,
        )
        wavefunction = Wavefunction(container)
        self.last_result = result
        self.last_energies = np.asarray([result.energy])
        return result.energy, wavefunction
