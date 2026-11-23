#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact thermal averages of a spin-1/2 XXZ chain, the reference of check_qmc_exact.py.

    H = sum_<ij> [Jz Sz_i Sz_j + Jxy/2 (S+_i S-_j + S-_i S+_j)] - h sum_i Sz_i

on L sites with periodic or open boundaries, by full diagonalization of the
2^L-dimensional Hamiltonian built here from the formula, not from ALPS's model
library. Prints the averages under the names the QMC applications use:
M = sum_i Sz_i, Ms = sum_i (-1)^i Sz_i, susceptibilities per site,
chi = (int_0^beta <M(tau) M> dtau - beta <M>^2) / L (the same for Ms), and the
specific heat per site. Needs NumPy; CTest does not run it.

Usage: exact_spin_chain.py L T Jz Jxy h periodic|open
"""
import sys

import numpy as np


def spin_operators(length):
    sz = np.diag([0.5, -0.5])
    sp = np.array([[0.0, 1.0], [0.0, 0.0]])

    def site(op, i):
        return np.kron(np.kron(np.eye(2**i), op), np.eye(2 ** (length - i - 1)))

    return [site(sz, i) for i in range(length)], [site(sp, i) for i in range(length)]


def averages(length, temperature, jz, jxy, h, periodic):
    sz, sp = spin_operators(length)
    bonds = [(i, (i + 1) % length) for i in range(length if periodic else length - 1)]
    ham = sum(jz * sz[i] @ sz[j] + 0.5 * jxy * (sp[i] @ sp[j].T + sp[i].T @ sp[j]) for i, j in bonds)
    ham = ham - h * sum(sz)
    beta = 1 / temperature
    w, v = np.linalg.eigh(ham)
    p = np.exp(-beta * (w - w.min()))
    p /= p.sum()

    def expectation(op):
        return float(p @ np.einsum("ij,ji->i", v.T @ op, v))

    def kubo(op):
        a = np.abs(v.T @ op @ v) ** 2
        dw = w[:, None] - w[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(np.abs(dw) < 1e-10, beta * p[:, None], (p[None, :] - p[:, None]) / dw)
        return float((a * f).sum())

    m = sum(sz)
    ms = sum((-1) ** i * s for i, s in enumerate(sz))
    energy = float(p @ w)
    return {
        "Energy Density": energy / length,
        "Specific Heat": beta**2 * (float(p @ w**2) - energy**2) / length,
        "Magnetization Density": expectation(m) / length,
        "Magnetization^2": expectation(m @ m),
        "Staggered Magnetization^2": expectation(ms @ ms),
        "Susceptibility": (kubo(m) - beta * expectation(m) ** 2) / length,
        "Staggered Susceptibility": (kubo(ms) - beta * expectation(ms) ** 2) / length,
    }


if __name__ == "__main__":
    L, T, Jz, Jxy, h, boundary = sys.argv[1:]
    for name, value in averages(int(L), float(T), float(Jz), float(Jxy), float(h), boundary == "periodic").items():
        print(f'"{name}": {value:.7g},')
