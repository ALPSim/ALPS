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
# docs/tutorials/time_evolution/domain_wall.py (itself a port of the legacy
# ALPS TEBD-02 "kink" tutorial), which evolved a product state by TEBD.
# Adaptation: rewritten for pyalps.classic_dmrg.evolve; the domain wall is
# prepared as the DMRG ground state of opposite fields on the two halves
# (with a weak exchange, see below); the exact references, chain length,
# time window and checks follow the alpsx script, the plots are condensed.

# %% [markdown]
# # Domain-wall ("kink") dynamics of an XX spin-1/2 chain
#
# An open chain starts in the domain wall
# $\lvert\uparrow\cdots\uparrow\,\downarrow\cdots\downarrow\rangle$ and
# evolves under the XX model
# $H = J_{xy}\sum_i (S^x_i S^x_{i+1} + S^y_i S^y_{i+1})$. The exchange melts
# the wall into a magnetization front that spreads ballistically.
#
# The XX chain maps onto free fermions, and the magnetization a distance
# $n$ right of the centre is a sum of Bessel functions [Antal *et al.*,
# Phys. Rev. E **59**, 4912 (1999)],
# $M(n,t) = -\tfrac12 \sum_{i=1-n}^{n-1} J_i(t)^2$, which at long times
# collapses onto $\phi(n/t) = -\tfrac1\pi \arcsin(n/t)$.

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv

from pyalps import classic_dmrg as dmrg
from pyalps import models


def exact_M(n, t):
    """M(n, t) = -1/2 sum_{i=1-n}^{n-1} J_i(t)^2."""
    return -0.5 * sum(jv(i, t) ** 2 for i in range(1 - n, n))


def scaling_phi(x):
    """-(1/pi) arcsin(min(x, 1))."""
    return -np.arcsin(min(x, 1.0)) / np.pi


L, DT, STEPS = 24, 0.05, 120  # t <= 6: the front stays away from the open ends
half = L // 2

# %% [markdown]
# ## Prepare the domain wall
#
# The wall is the ground state of fields $\mp S^z_i$ on the left and right
# halves. Its DMRG density matrices would have rank one, and the warm-up of
# the classic DMRG can then lose the charge sectors it needs, so a weak
# exchange ($10^{-3}$) is added: it changes the local magnetization by
# about $10^{-7}$ and lets the sweeps find the wall.

# %%
lattice = models.lattice.chain(L, bc="open")
model = models.xxz(lattice, Jxy=1.0, Jz=0.0)
site = model.site_for(0)
wall = models.TermSum()
for i in range(L):
    wall += (-1.0 if i < half else 1.0) * site.Sz(i)
for i in range(L - 1):
    wall += 0.5e-3 * (site.Sp(i) * site.Sm(i + 1) + site.Sm(i) * site.Sp(i + 1))

# %% [markdown]
# ## Evolve under the XX chain
#
# `measurements={"Sz": "Sz"}` measures $S^z$ on every site; `measure_steps=4`
# measures every fourth step. The kept states follow a truncation error of
# $10^{-10}$, up to 40 states.

# %%
run = dmrg.evolve(wall, model, time_step=DT, steps=STEPS, states=40, truncation_error=1e-10,
                  sector={"Sz": 0}, measure_steps=4, measurements={"Sz": "Sz"})
ts = run.times
profiles = run.measurements["Sz"].values
print(f"measured {len(ts)} times up to t = {ts[-1]:.1f}; discarded weight {run.truncation_errors.sum():.1e}")
print("t = 0 profile:", np.round(profiles[0], 3))

# %% [markdown]
# ## Against the exact solution

# %%
print(f"{'t':>6} {'n':>2} {'M_numeric':>11} {'M_exact':>11} {'|diff|':>9}")
max_err = {1: 0.0, 2: 0.0}
for k, t in enumerate(ts):
    for n in (1, 2):
        m_num, m_ex = profiles[k, half + n - 1], exact_M(n, t)
        if t <= 4.0:
            max_err[n] = max(max_err[n], abs(m_num - m_ex))
        if k % 6 == 0:
            print(f"{t:6.2f} {n:>2} {m_num:11.5f} {m_ex:11.5f} {abs(m_num - m_ex):9.2e}")

scaling = [(n / t, profiles[k, half + n - 1], scaling_phi(n / t))
           for k, t in enumerate(ts) if t >= 2.0 for n in range(1, 6) if n / t <= 1.2]
x, numeric, phi = (np.array(column) for column in zip(*scaling))
rmse = float(np.sqrt(np.mean((numeric - phi) ** 2)))

front = np.array([np.max(np.nonzero(p[half:] > -0.4)[0], initial=-1) + 1 for p in profiles])

# %% [markdown]
# ## Checks
#
# The initial state is the wall and the total $S^z$ stays zero; the centre
# melts while the ends are untouched and the front advances; the profile
# follows the Bessel solution at $n = 1, 2$ (the remaining error is the
# Trotter error, $O(dt^2)$) and approaches the scaling form.

# %%
assert np.allclose(profiles[0, :half], 0.5, atol=1e-4) and np.allclose(profiles[0, half:], -0.5, atol=1e-4)
assert np.all(np.abs(profiles.sum(axis=1)) < 1e-8), "total Sz must stay zero"
assert abs(profiles[-1, half]) < abs(profiles[0, half]) - 0.2, "the centre must melt"
assert abs(profiles[-1, 0]) > 0.45 and abs(profiles[-1, -1]) > 0.45, "the ends must be untouched"
assert front[-1] > front[1], "the front must advance"
assert max_err[1] < 1e-3 and max_err[2] < 1e-3, max_err
assert rmse < 0.08, rmse
print(f"max |M_num - M_exact| (t <= 4): n=1 {max_err[1]:.1e}, n=2 {max_err[2]:.1e}; scaling RMSE {rmse:.3f}")
print("All physics assertions passed.")

# %%
fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(14.5, 4.3))
distance = np.arange(L) - (half - 0.5)
for k in np.linspace(0, len(ts) - 1, 6).astype(int):
    ax_a.plot(distance, profiles[k], label=f"$t={ts[k]:.1f}$")
ax_a.set(xlabel="distance from the centre", ylabel=r"$\langle S^z_i\rangle$", title="Front spreading")
ax_a.legend(fontsize=7, ncol=2)
for n, color in ((1, "C0"), (2, "C3")):
    ax_b.plot(ts, profiles[:, half + n - 1], "o", ms=4, color=color, label=f"tDMRG $n={n}$")
    ax_b.plot(ts, [exact_M(n, t) for t in ts], "-", color=color, alpha=0.7, label=f"exact $n={n}$")
ax_b.set(xlabel="time $t$", ylabel="$M(n,t)$", title="Against the Bessel solution")
ax_b.legend(fontsize=8)
order = np.argsort(x)
ax_c.plot(x[order], numeric[order], "o", ms=4, color="C2", label="tDMRG")
grid = np.linspace(0, 1.2, 100)
ax_c.plot(grid, [scaling_phi(v) for v in grid], "k-", label=r"$-\frac{1}{\pi}\arcsin(n/t)$")
ax_c.set(xlabel="$n/t$", ylabel="$M(n,t)$", title="Scaling collapse", xlim=(0, 1.2))
ax_c.legend(fontsize=8)
fig.tight_layout()
plt.show()
