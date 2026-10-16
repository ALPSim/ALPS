# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/src/qdk_chemistry_alps/ed_calculator.py. Changes:
# the eigenstates come from pyalps.ed through abinitio.solve_ground_state, and
# the "engine" setting, whose values all selected the one native solver, is
# gone.
"""AlpsEdCalculator: qdk-chemistry MultiConfigurationCalculator backed by
pyalps's exact diagonalization of the operator-algebra Hamiltonian.

Input : qdk_chemistry.data.Hamiltonian wrapping supported real, restricted
        one-/two-body integrals in chemist notation + (n_alpha, n_beta).
Output: (total energy, Wavefunction wrapping a StateVectorContainer with the
        CI vector and determinant list). The default keeps all amplitudes;
        a positive ci_coefficient_cutoff drops small ones and renormalizes.

Determinant conventions: ALPS basis rows are per-site local states with the
interleaved dn-before-up mode ordering; qdk Configurations use per-orbital
'0/u/d/2' strings (leftmost char = orbital 0) with the blocked
alpha-block-then-beta-block creation ordering (PySCF-style). The parity map
in abinitio.py converts coefficients between the two.
"""

from __future__ import annotations

import numpy as np
from qdk_chemistry.algorithms import MultiConfigurationCalculator
from qdk_chemistry.data import (
    Configuration,
    Settings,
    StateVectorContainer,
    Wavefunction,
)

from .abinitio import (
    alps_state_to_occupations,
    build_abinitio_hamiltonian,
    interleaved_to_blocked_parity,
    solve_ground_state,
)

__all__ = ["AlpsEdCalculator"]

_STATE_CHAR = {(0, 0): "0", (1, 0): "u", (0, 1): "d", (1, 1): "2"}


class AlpsEdCalculator(MultiConfigurationCalculator):
    """Exact diagonalization via pyalps's operator algebra and pyalps.ed."""

    def __init__(self):
        super().__init__()
        self._settings = Settings()
        self._settings._set_default("integral_cutoff", "double", 1e-12)
        self._settings._set_default("ci_coefficient_cutoff", "double", 0.0)

    def name(self) -> str:
        return "alps_ed"

    def _run_impl(self, hamiltonian, n_alpha: int, n_beta: int):
        h1a, h1b = hamiltonian.get_one_body_integrals()
        h1a = np.asarray(h1a)
        if not hamiltonian.is_restricted() and not np.allclose(h1a, np.asarray(h1b)):
            raise NotImplementedError(
                "alps_ed currently supports restricted (spin-symmetric) Hamiltonians only"
            )
        n = h1a.shape[0]

        two = hamiltonian.get_two_body_integrals()
        eri = np.asarray(two[0] if isinstance(two, (tuple, list)) else two)
        eri = eri.reshape(n, n, n, n)
        ecore = float(hamiltonian.get_core_energy())

        H = build_abinitio_hamiltonian(
            h1a, eri, tol=float(self._settings.get("integral_cutoff"))
        )
        res = solve_ground_state(H, n_alpha, n_beta)
        energy = float(res.eigenvalues[0]) + ecore
        vec = np.asarray(res.eigenvectors)[:, 0]
        states = np.asarray(res.problem.states)

        cutoff = float(self._settings.get("ci_coefficient_cutoff"))
        coeffs, dets = [], []
        order = np.argsort(-np.abs(vec))
        for i in order:
            c = float(vec[i])
            if cutoff > 0.0 and abs(c) < cutoff:
                break
            a_occ, b_occ = alps_state_to_occupations(states[i])
            c *= interleaved_to_blocked_parity(a_occ, b_occ)
            s = "".join(_STATE_CHAR[(int(a), int(b))] for a, b in zip(a_occ, b_occ, strict=True))
            dets.append(Configuration.from_spin_half_string(s))
            coeffs.append(c)

        coeffs = np.asarray(coeffs)
        coeffs /= np.linalg.norm(coeffs)
        container = StateVectorContainer(coeffs, dets, hamiltonian.get_orbitals())
        return energy, Wavefunction(container)
