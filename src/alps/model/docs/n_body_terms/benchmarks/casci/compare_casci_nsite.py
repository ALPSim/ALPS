"""Molecular Hamiltonians, full (pq|rs), QDK CASCI vs ALPS ED with n-site terms.

    ALPS_BIN=<alps>/_build/default/bin ALPS_XML_PATH=<alps>/src/alps/resources \
    python compare_casci_nsite.py        # needs qdk-chemistry

compare_ed.py covers lattice models, whose two-body part fits SITETERM and
BONDTERM. A molecular Hamiltonian does not: its (pq|rs) touches 3 and 4
distinct orbitals. qdk2alps now emits those as <PLAQUETTE>s plus
<PLAQUETTETERM>s, with nothing dropped, so ALPS ED must reproduce CASCI.

Each case: SCF -> MO integrals (all orbitals active) -> one Hamiltonian ->
QDK macis_cas and ALPS fulldiag. A third column diagonalises the same
integrals with a dense Jordan-Wigner matrix in numpy, independent of both.
"""

import itertools
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import qdk_alps_ed                                    # noqa: E402,F401  registers alps_ed
from qdk2alps import GENERAL_MODEL, former_for, plaquette_groups  # noqa: E402
from qdk_chemistry.algorithms import create            # noqa: E402
from qdk_chemistry.data import Element, Structure      # noqa: E402

TOL = 1e-8
ANG = 1.8897259886                                     # Bohr per Angstrom


def molecule(coords_ang, elements, charge=0, mult=1, basis="sto-3g"):
    s = Structure(np.asarray(coords_ang, float) * ANG, elements)
    with open(os.devnull, "w") as null:                # SCF prints a report
        fd = os.dup(1)
        os.dup2(null.fileno(), 1)
        try:
            _, wf = create("scf_solver").run(s, charge, mult, basis)
            ham = create("hamiltonian_constructor").run(wf.get_orbitals())
        finally:
            os.dup2(fd, 1)
            os.close(fd)
    return ham


def dense_reference(ham, na, nb):
    """Sector spectrum from a dense Jordan-Wigner matrix. Modes 2p+spin."""
    h1 = np.asarray(ham.get_one_body_integrals()[0])
    n = h1.shape[0]
    g = np.asarray(ham.get_two_body_integrals()[0]).reshape(n, n, n, n)
    M = 2 * n
    states = [s for s in range(1 << M)
              if bin(s & int("01" * n, 2)).count("1") == na
              and bin(s & int("10" * n, 2)).count("1") == nb]
    index = {s: i for i, s in enumerate(states)}

    def apply(ops, s):           # ops applied right to left: list of (dag, mode)
        sign = 1
        for dag, m in reversed(ops):
            occ = (s >> m) & 1
            if occ == dag:
                return None, 0
            if bin(s & ((1 << m) - 1)).count("1") % 2:
                sign = -sign
            s ^= 1 << m
        return s, sign

    H = np.zeros((len(states), len(states)))
    for col, s in enumerate(states):
        for p, q in itertools.product(range(n), repeat=2):
            if h1[p, q] == 0:
                continue
            for a in (0, 1):
                t, sg = apply([(1, 2 * p + a), (0, 2 * q + a)], s)
                if t is not None:
                    H[index[t], col] += h1[p, q] * sg
        for p, q, r, s_ in zip(*np.nonzero(np.abs(g) > 1e-14)):
            v = 0.5 * g[p, q, r, s_]
            for a in (0, 1):
                for b in (0, 1):
                    t, sg = apply([(1, 2 * p + a), (1, 2 * r + b),
                                   (0, 2 * s_ + b), (0, 2 * q + a)], s)
                    if t is not None:
                        H[index[t], col] += v * sg
    return np.linalg.eigvalsh(H) + float(ham.get_core_energy())


def cases():
    H = Element.H
    yield "H4 chain 0.9 A", molecule([[0, 0, 0.9 * i] for i in range(4)], [H] * 4), [(2, 2), (2, 1), (3, 1)]
    yield "H4 rectangle 0.8x1.2 A", molecule([[0, 0, 0], [0.8, 0, 0], [0.8, 1.2, 0], [0, 1.2, 0]], [H] * 4), [(2, 2), (3, 1)]
    yield "H4 bent/asym", molecule([[0, 0, 0], [0.75, 0.1, 0], [1.4, 0.9, 0.2], [2.4, 1.0, -0.3]], [H] * 4), [(2, 2), (2, 1), (1, 1)]
    yield "H3+ triangle", molecule([[0, 0, 0], [0.9, 0, 0], [0.45, 0.8, 0]], [H] * 3, charge=1), [(1, 1), (2, 0)]
    yield "H4 chain 1.5 A (stretched)", molecule([[0, 0, 1.5 * i] for i in range(4)], [H] * 4), [(2, 2), (3, 1)]
    yield "LiH 1.6 A (6 orbitals)", molecule([[0, 0, 0], [0, 0, 1.6]], [Element.Li, H]), [(2, 2), (2, 1)]
    yield "H6 ring 1.0 A", molecule([[np.cos(k * np.pi / 3), np.sin(k * np.pi / 3), 0] for k in range(6)], [H] * 6), [(3, 3)]


def main():
    casci = create("multi_configuration_calculator")
    alps = create("multi_configuration_calculator", "alps_ed")
    print(f"{'case':30}{'sector':>8}{'groups 2/3/4':>14}{'QDK CASCI':>16}{'ALPS ED':>16}"
          f"{'|dE|':>9}{'dense JW':>16}{'spec |dE|':>10}{'t_alps':>8}")
    print("-" * 127)
    bad = 0
    for label, ham, sectors in cases():
        h1 = np.asarray(ham.get_one_body_integrals()[0])
        model, _ = former_for(ham, h1)
        groups = plaquette_groups(ham, h1.shape[0])
        sizes = [sum(len(k) == m for k in groups) for m in (2, 3, 4)]
        for na, nb in sectors:
            e_qdk, _ = casci.run(ham, na, nb)
            t0 = time.time()
            e_alps, _ = alps.run(ham, na, nb)
            dt = time.time() - t0
            levels = np.asarray(alps.last_levels)
            ref = dense_reference(ham, na, nb)
            dspec = np.max(np.abs(np.sort(levels) - ref)) if len(levels) == len(ref) else np.inf
            ok = model == GENERAL_MODEL and abs(e_qdk - e_alps) < TOL and dspec < TOL
            bad += not ok
            print(f"{label:30}{f'({na},{nb})':>8}{'/'.join(map(str, sizes)):>14}"
                  f"{e_qdk:16.10f}{e_alps:16.10f}{abs(e_qdk - e_alps):9.1e}"
                  f"{ref[0]:16.10f}{dspec:10.1e}{dt:7.1f}s" + ("" if ok else "  <-- MISMATCH"))
    print("-" * 127)
    print("ALL AGREE" if not bad else f"{bad} MISMATCH(ES)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
