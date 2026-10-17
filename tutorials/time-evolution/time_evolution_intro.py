# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# Copyright (C) 2026 ALPS Collaboration
# SPDX-License-Identifier: MIT
#
# Ported from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227, file
# docs/tutorials/time_evolution/time_evolution_intro.py, which ran TEBD on
# matrix product states (alps.mps.evolve, imaginary then real time from a
# random state). Adaptation: rewritten for pyalps.classic_dmrg.evolve, the
# adaptive time-dependent DMRG of the classic engine. Every evolution here
# starts from the DMRG ground state of an initial Hamiltonian, so the two
# regimes are two runs from the XX chain's ground state; the Trotter order
# is two, and the energy-drift check compares two time steps.

# %% [markdown]
# # Time evolution basics (imaginary and real time)
#
# `pyalps.classic_dmrg.evolve` evolves a state by **adaptive time-dependent
# DMRG** (White and Feiguin, PRL 93, 076401 (2004)): the evolution operator
# is split into second-order Trotter gates on the bonds,
# $e^{-dt H} \approx e^{-dt A/2} e^{-dt B} e^{-dt A/2}$ with $A$ the even
# bonds and $B$ the odd ones, and each gate is applied while the DMRG sweeps
# over its bond, so the kept basis follows the state.
#
# The initial state is always the DMRG ground state of an *initial*
# Hamiltonian (in a charge sector); it evolves under a *final* Hamiltonian
# on the same sites, which the Trotter gates restrict to neighbouring sites
# in chain order. (`method="targeting"` lifts that restriction and the
# Trotter error, by time-step targeting with a Krylov step.)
# There are two regimes:
#
# * **imaginary time** (`imaginary=True`) applies $e^{-\tau H}$ and
#   renormalizes: every excited component decays relative to the lowest
#   state the initial state overlaps;
# * **real time** applies $e^{-iHt}$: quenches, transport, echoes.

# %%
import numpy as np
import matplotlib.pyplot as plt

from pyalps import classic_dmrg as dmrg
from pyalps import models

lattice = models.lattice.chain(8, bc="open")
xx = models.xxz(lattice, Jxy=1.0, Jz=0.0)       # the initial Hamiltonian
heisenberg = models.heisenberg(lattice, J=1.0)  # the final one
SECTOR = {"Sz": 0}
EXACT = {"sweeps": dmrg.fixed_sweeps(4, 64), "warmup_states": 64}  # exact on eight sites

# %% [markdown]
# ## Imaginary time: projecting onto the Heisenberg ground state
#
# Starting from the XX chain's ground state, $e^{-\tau H}$ with the
# Heisenberg Hamiltonian converges to the Heisenberg ground state. The
# energy decreases monotonically and settles $O(dt^2)$ above the exact
# value, at the fixed point of the Trotterized propagator.

# %%
e0 = dmrg.solve(heisenberg, sector=SECTOR, **EXACT).energy
cooling = dmrg.evolve(xx, heisenberg, time_step=0.05, steps=100, imaginary=True, states=64,
                      measure_steps=5, sector=SECTOR, **EXACT)
print(f"DMRG ground-state energy          {e0:+.8f}")
print(f"after imaginary time tau = {cooling.times[-1]:.1f}   {cooling.energies[-1]:+.8f}")

# %% [markdown]
# ## Real time: a quench from the XX chain to the Heisenberg chain
#
# The XX ground state is not an eigenstate of the Heisenberg chain, so it
# evolves. The energy $\langle H \rangle$ is conserved up to the Trotter
# error, which is second order in the step: halving the step quarters the
# drift.

# %%
drifts = {}
for dt in (0.1, 0.05):
    quench = dmrg.evolve(xx, heisenberg, time_step=dt, steps=round(5.0 / dt), states=64,
                         measure_steps=round(0.25 / dt), sector=SECTOR, **EXACT)
    drifts[dt] = (quench, np.abs(quench.energies - quench.energies[0]).max())
    print(f"dt = {dt:<5} max |E(t) - E(0)| = {drifts[dt][1]:.3e}")
quench = drifts[0.05][0]

# %% [markdown]
# ## Checks

# %%
assert np.all(np.diff(cooling.energies) < 1e-8), "imaginary time must lower the energy"
assert abs(cooling.energies[-1] - e0) < 5e-3, "imaginary time must reach the ground state"
assert cooling.energies[0] > e0 + 0.1, "the XX ground state is not the Heisenberg one"
ratio = drifts[0.1][1] / drifts[0.05][1]
assert 2.5 < ratio < 6.0, f"the energy drift should shrink as dt^2 (ratio {ratio:.2f})"
assert drifts[0.05][1] < 1e-3
print("All physics assertions passed.")

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
axes[0].plot(cooling.times, cooling.energies, "o-")
axes[0].axhline(e0, color="0.5", ls=":", label="DMRG ground state")
axes[0].set(xlabel=r"imaginary time $\tau$", ylabel="energy", title="Imaginary time")
axes[0].legend()
for dt, (run, _) in drifts.items():
    axes[1].plot(run.times, run.energies - run.energies[0], "o-", ms=3, label=f"dt = {dt}")
axes[1].set(xlabel="time $t$", ylabel=r"$E(t) - E(0)$", title="Real time")
axes[1].legend()
for ax in axes:
    ax.grid(alpha=0.3)
fig.tight_layout()
plt.show()
