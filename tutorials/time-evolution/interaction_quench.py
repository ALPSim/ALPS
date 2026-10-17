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
# docs/tutorials/time_evolution/interaction_quench.py (itself a port of the
# legacy ALPS TEBD-01 "bhquench" tutorial 1a-1c), which ran TEBD from a
# DMRG state with an additive drive sampled at step endpoints. Adaptation:
# rewritten for pyalps.classic_dmrg.evolve, which prepares the initial state
# itself and samples a drive at step midpoints; the echo comes from the
# kept states (loschmidt_echo); protocols, parameters and checks follow the
# alpsx script.

# %% [markdown]
# # Interaction quench of a hard-core boson chain (Loschmidt echo)
#
# The open chain of hard-core bosons ($L = 10$, $N = 5$)
#
# $$
# H = -t \sum_i \bigl(b^\dagger_i b_{i+1} + b_i b^\dagger_{i+1}\bigr) + V \sum_i n_i n_{i+1}
# $$
#
# is a charge-density wave for $V/t \gg 1$ and a superfluid for small
# $V/t$. Starting from the $V = 10$ ground state we drive $V(t)$ down and
# back up and watch the **Loschmidt echo**
# $\mathcal{L}(t) = |\langle\psi(0)|\psi(t)\rangle|^2$; how close it returns
# to one measures how adiabatic the round trip was.
#
# * **1a**, linear ramp $10 \to 0 \to 10$, fast against slow;
# * **1b**, ramp down, hold at $V = 0$, ramp up, without and with a hold;
# * **1c**, power-law ramps $V \sim (t/\tau)^p$ with $p = 1$ and $p = 3$.
#
# The final Hamiltonian is the chain at $V = 0$; the interaction enters as
# a *drive*: `drives=[(interaction, V)]` adds $V(t) \sum_i n_i n_{i+1}$, with
# the function $V(t)$ sampled at the middle of every step.

# %%
import numpy as np
import matplotlib.pyplot as plt

from pyalps import classic_dmrg as dmrg
from pyalps import models

L, HALF, VMAX, DT = 10, 5, 10.0, 0.1
lattice = models.lattice.chain(L, bc="open")
initial = models.hardcore_boson_model(lattice, t=1.0, V=VMAX)
free = models.hardcore_boson_model(lattice, t=1.0, V=0.0)
interaction = models.hardcore_boson_model(lattice, t=0.0, V=1.0)  # sum_i n_i n_{i+1}


def ramp_linear(tau):
    """1a: V 10 -> 0 over [0, tau], then 0 -> 10 over [tau, 2 tau]."""
    return lambda t: VMAX * (1.0 - t / tau) if t < tau else VMAX * min((t - tau) / tau, 1.0)


def ramp_hold(tau_ramp, tau_hold):
    """1b: 10 -> 0 over tau_ramp, hold 0 for tau_hold, 0 -> 10 over tau_ramp."""
    def f(t):
        if t < tau_ramp:
            return VMAX * (1.0 - t / tau_ramp)
        if t < tau_ramp + tau_hold:
            return 0.0
        return VMAX * min((t - tau_ramp - tau_hold) / tau_ramp, 1.0)
    return f


def ramp_power(tau, p):
    """1c: V 10 -> 0 -> 10, each half as (t / tau)^p."""
    return lambda t: VMAX * (1.0 - (t / tau) ** p) if t < tau else VMAX * min((t - tau) / tau, 1.0) ** p


def quench(ramp, total_t):
    """One round trip from the V = 10 ground state: (t, echo, V(t))."""
    steps = max(2, round(total_t / DT))
    run = dmrg.evolve(initial, free, time_step=DT, steps=steps, drives=[(interaction, ramp)], states=40,
                      truncation_error=1e-9, measure_steps=2, keep_mps=True, sector={"N": HALF})
    return run.times, dmrg.loschmidt_echo(run), np.array([ramp(t) for t in run.times])


# %% [markdown]
# ## The three protocols

# %%
runs = {"1a": {}, "1b": {}, "1c": {}}
for tau, label in ((0.5, "fast (tau=0.5)"), (4.0, "slow (tau=4.0)")):
    runs["1a"][label] = quench(ramp_linear(tau), 2.0 * tau)
for hold, label in ((0.0, "no hold"), (4.0, "hold=4.0")):
    runs["1b"][label] = quench(ramp_hold(2.0, hold), 4.0 + hold)
for p, label in ((1.0, "p=1.0 (linear)"), (3.0, "p=3.0 (convex)")):
    runs["1c"][label] = quench(ramp_power(3.0, p), 6.0)
for name, sub in runs.items():
    for label, (_, echo, _) in sub.items():
        print(f"{name} {label:>16}: echo opens {echo[0]:.3f}, min {echo.min():.3f}, ends {echo[-1]:.3f}")

# %% [markdown]
# ## Checks
#
# Every echo opens at one and stays in $[0, 1]$; the slow round trip of 1a
# returns closer to the initial state than the fast one; a hold at $V = 0$
# does not raise the echo's minimum; the two ramp shapes of 1c differ.

# %%
for sub in runs.values():
    for label, (_, echo, _) in sub.items():
        assert abs(echo[0] - 1.0) < 1e-10, label
        assert echo.min() >= -1e-9 and echo.max() <= 1.0 + 1e-6, label
fast_end, slow_end = runs["1a"]["fast (tau=0.5)"][1][-1], runs["1a"]["slow (tau=4.0)"][1][-1]
assert slow_end > fast_end + 0.05, (slow_end, fast_end)
assert runs["1b"]["hold=4.0"][1].min() <= runs["1b"]["no hold"][1].min() + 1e-3
linear, convex = runs["1c"]["p=1.0 (linear)"][1], runs["1c"]["p=3.0 (convex)"][1]
assert np.max(np.abs(linear - convex)) > 0.02
print("All physics assertions passed.")

# %%
fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.6), gridspec_kw={"height_ratios": [2.0, 1.0]})
titles = {"1a": r"1a  linear ramp (sweep $\tau$)", "1b": "1b  hold at $V=0$", "1c": "1c  power-law ramp"}
for column, name in enumerate(runs):
    for (label, (t, echo, v)), color in zip(runs[name].items(), ("C0", "C3")):
        axes[0][column].plot(t, echo, "o-", ms=3, color=color, label=label)
        axes[1][column].plot(t, v, color=color)
    axes[0][column].set(title=titles[name], ylim=(-0.02, 1.02))
    axes[0][column].legend(fontsize=8)
    axes[1][column].set(xlabel="time $t$")
axes[0][0].set_ylabel(r"$|\langle\psi(0)|\psi(t)\rangle|^2$")
axes[1][0].set_ylabel("$V(t)$")
fig.tight_layout()
plt.show()
