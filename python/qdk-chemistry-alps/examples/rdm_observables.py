# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/rdm_observables.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 483100f45061e1af3d36233dd5461b5d87553790fbceb23e24fe6440ada1c4cc.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Fourteen-site Hubbard observables without constructing a CI vector.

This is explicitly an ALPS-native backup. Returning a QDK MPSContainer is
demonstrated separately by large_active_space.py with the optional QDK build.
"""

import time
from math import comb

from pyalps import models
import numpy as np
from pyalps import classic_dmrg as dmrg
from pyalps.models import hubbard
from demo_support import finish, parser


def plot(plt, data):
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    axes[0].plot(data["natural_occupations"], "o-", color="#087e8b")
    axes[0].set(xlabel="Natural orbital (descending occupation)", ylabel="Occupation", ylim=(0, 2),
                title="14-site Hubbard: one-particle RDM")
    axes[1].plot(range(14), data["double_occupancy"], "o-", color="#d68c22")
    axes[1].set(xlabel="Site", ylabel="⟨n↑ n↓⟩", title="Local two-particle information")
    return figure


def main():
    args = parser(__doc__).parse_args()
    started = time.perf_counter()
    h = hubbard(models.lattice.chain(14, bc="open"), t=1., U=4.)
    from qdk_chemistry_alps.interop import hamiltonian_to_integrals
    h1, eri, core = hamiltonian_to_integrals(h)
    result = dmrg.solve_molecular(h1, eri, 7, 7, core_energy=core,
                                 sweeps=dmrg.fixed_sweeps(8,64), lanczos_tolerance=1e-10)
    gamma = result.one_rdm
    doubles = dmrg.local_expectations(result, "n_updn")
    # For Hubbard only these two-particle diagonal elements are needed.
    from_rdms = -np.trace(gamma, offset=1) - np.trace(gamma, offset=-1) + 4 * doubles.sum()
    np.testing.assert_allclose(from_rdms, result.energy, atol=1e-7, rtol=0)
    np.testing.assert_allclose(np.trace(gamma), 14, atol=1e-8, rtol=0)
    natural = np.linalg.eigvalsh(gamma)[::-1]
    assert np.min(natural) > -1e-9 and np.max(natural) < 2 + 1e-9
    determinants = comb(14, 7)**2
    print("ALPS-native MPS result; no returned QDK RdmContainer and no CI reconstruction")
    print(f"14-site Hubbard: E/t={result.energy:.10f}; RDM energy error={abs(from_rdms-result.energy):.2e}")
    print(f"{determinants:,} determinants would need {8*determinants/1e6:.1f} MB for real amplitudes alone.")
    print("This checks energy consistency and particle count, not DMRG convergence.")
    finish(args, {"model": {"sites": 14, "t": 1., "U": 4., "sector": {"N": 14, "Sz": 0}, "boundary": "open"},
                  "controls": {"bond": 64, "sweeps": 8,  "truncation": 0.0, "eigensolver_tol": 1e-10},
                  "workflow": "ALPS-native observables; not a QDK RdmContainer",
                  "energy": result.energy, "energy_from_rdms": float(from_rdms.real),
                  "natural_occupations": natural, "double_occupancy": doubles,
                  "one_rdm": gamma.real, "determinants": determinants}, started, plot)


if __name__ == "__main__":
    main()
