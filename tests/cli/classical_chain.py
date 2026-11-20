# SPDX-License-Identifier: MIT
"""Exact thermal averages of an open chain of classical unit spins.

For H = -J sum_i s_i.s_{i+1} on an open chain the relative orientations of
neighbouring spins are independent, so with K = J/T every bond has the same
<c> = <s_i.s_{i+1}> and <c^2>, and <s_i.s_j> = <c>^|i-j|:

    Ising       <c> = tanh K,            <c^2> = 1
    XY          <c> = I1(K)/I0(K),       <c^2> = (1 + I2(K)/I0(K)) / 2
    Heisenberg  <c> = coth K - 1/K,      <c^2> = 1 - 2 <c> / K
"""
import math


def _bessel_i(n, k):
    """Modified Bessel function I_n(k) from its power series."""
    return sum((k / 2) ** (2 * m + n) / (math.factorial(m) * math.factorial(m + n)) for m in range(60))


def bond(model, k):
    """<c> and <c^2> of one bond at K = J/T."""
    if model == "Ising":
        return math.tanh(k), 1.0
    if model == "XY":
        i0 = _bessel_i(0, k)
        return _bessel_i(1, k) / i0, (1 + _bessel_i(2, k) / i0) / 2
    if model == "Heisenberg":
        c = 1 / math.tanh(k) - 1 / k
        return c, 1 - 2 * c / k
    raise ValueError(model)


def magnetization_squared(model, k, length, staggered=False):
    """<M^2> / L^2 with M = sum_i s_i, or sum_i (-1)^i s_i if staggered."""
    c = bond(model, k)[0] * (-1 if staggered else 1)
    return sum(c ** abs(i - j) for i in range(length) for j in range(length)) / length**2


def energy_density(model, k, coupling, length):
    """<H> / L."""
    return -coupling * (length - 1) * bond(model, k)[0] / length


def specific_heat(model, k, length):
    """(<H^2> - <H>^2) / (T^2 L); the bond terms are independent."""
    c, c2 = bond(model, k)
    return k**2 * (length - 1) * (c2 - c**2) / length
