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
# docs/tutorials/time_evolution/quench_hardcore_boson.py, which prepared the
# CDW with DMRG and evolved it by fourth-order TEBD. Adaptation: rewritten
# for pyalps.classic_dmrg.evolve (second-order Trotter, the state prepared
# by the call itself, drives sampled at step midpoints, the echo from the
# kept states); the figures are condensed and the ramped part gains checks.

# %% [markdown]
# # Sudden and ramped quenches of the hard-core boson chain
#
# The hard-core boson chain
# $H = -t \sum_i (b^\dagger_i b_{i+1} + \mathrm{h.c.}) + V \sum_i n_i n_{i+1}$
# at half filling is a charge-density wave (CDW) for $V/t \gg 1$ and a
# superfluid for $V/t \ll 1$. We prepare the CDW at $V/t = 10$, switch the
# repulsion off and watch the density pattern melt; then we ramp $V$ down
# instead and compare fast and slow ramps through the Loschmidt echo.

# %%
import numpy as np
import matplotlib.pyplot as plt

from pyalps import classic_dmrg as dmrg
from pyalps import models

L, HALF = 10, 5
lattice = models.lattice.chain(L, bc="open")
cdw = models.hardcore_boson_model(lattice, t=1.0, V=10.0)
free = models.hardcore_boson_model(lattice, t=1.0, V=0.0)
SECTOR = {"N": HALF}

# %% [markdown]
# ## A sudden quench $V/t = 10 \to 0$
#
# `evolve` solves the CDW ground state (its result is `run.initial`) and
# evolves it under the free chain. Five bosons on ten open sites avoid every
# neighbour pair in C(6, 5) = 6 ways: one domain wall between the two CDW
# orders, at any of six places (the two perfect CDWs put it at an end).
# Hopping moves the wall, so the ground state spreads it over the chain:
# the density alternates strongly at the ends and is 1/2 in the middle.
# The total density is conserved by the dynamics, a check on the evolution.

# %%
run = dmrg.evolve(cdw, free, time_step=0.02, steps=500, states=40, truncation_error=1e-10,
                  measure_steps=10, measurements={"n": "n"}, sector=SECTOR)
times, density = run.times, run.measurements["n"].values
print(f"CDW ground-state energy {run.initial.energy:+.6f}")
print("initial density", np.round(density[0], 3))
print(f"evolved to t = {times[-1]:.1f}; N in [{density.sum(axis=1).min():.8f}, {density.sum(axis=1).max():.8f}]")
print("final density  ", np.round(density[-1], 3))

# %% [markdown]
# ## A ramped quench and the Loschmidt echo
#
# Now $V$ is ramped down linearly over a time $\tau$ and then held at zero;
# both ramps run for the same total time. A near-sudden ramp leaves the
# state at once and the echo $|\langle\psi(0)|\psi(t)\rangle|^2$ collapses;
# a slow ramp lets the state follow the instantaneous ground state for
# longer, so the echo stays high until $V$ has dropped.

# %%
interaction = models.hardcore_boson_model(lattice, t=0.0, V=1.0)  # sum_i n_i n_{i+1}


def ramp_down(tau):
    return lambda t: 10.0 * max(0.0, 1.0 - t / tau)


total_t, steps = 4.0, 60
echoes = {}
for tau in (0.5, 4.0):
    ramped = dmrg.evolve(cdw, free, time_step=total_t / steps, steps=steps, drives=[(interaction, ramp_down(tau))],
                         states=40, truncation_error=1e-9, measure_steps=2, keep_mps=True, sector=SECTOR)
    echoes[tau] = (ramped.times, dmrg.loschmidt_echo(ramped))
    print(f"tau = {tau}: echo opens {echoes[tau][1][0]:.4f}, min {echoes[tau][1].min():.4f}, "
          f"ends {echoes[tau][1][-1]:.4f}")

# %% [markdown]
# ## Checks

# %%
np.testing.assert_allclose(density.sum(axis=1), HALF, atol=1e-8)
# The CDW's domain wall is spread over the chain (above): the density
# alternates strongly towards both ends.
assert density[0, 0] > 0.85 and density[0, 1] < 0.15 and density[0, -1] > 0.85 and density[0, -2] < 0.15
contrast = lambda n: np.abs(np.diff(n)).mean()  # noqa: E731
assert contrast(density[-1]) < 0.5 * contrast(density[0]), "the CDW must melt"
(t_fast, fast), (t_slow, slow) = echoes[0.5], echoes[4.0]
assert abs(fast[0] - 1.0) < 1e-10 and abs(slow[0] - 1.0) < 1e-10
middle = np.searchsorted(t_fast, 1.0)  # after the fast ramp, during the slow one
assert slow[middle] > fast[middle] + 0.2, "the slow ramp must keep the echo higher"
print("All physics assertions passed.")

# %%
fig, (ax_n, ax_e) = plt.subplots(1, 2, figsize=(12, 4.4))
image = ax_n.imshow(density.T, aspect="auto", origin="lower", extent=[times[0], times[-1], 0.5, L + 0.5],
                    cmap="viridis", vmin=0.0, vmax=1.0)
ax_n.set(xlabel="time $t$", ylabel="site", title=r"Sudden quench $V/t = 10 \to 0$")
fig.colorbar(image, ax=ax_n, label=r"$\langle n_i\rangle$")
for tau, (t, echo) in echoes.items():
    ax_e.plot(t, echo, "o-", ms=3, label=rf"$\tau = {tau:g}$")
ax_e.set(xlabel="time $t$", ylabel=r"$|\langle\psi(0)|\psi(t)\rangle|^2$", title="Ramped quench", ylim=(-0.02, 1.02))
ax_e.legend()
ax_e.grid(alpha=0.3)
fig.tight_layout()
plt.show()
