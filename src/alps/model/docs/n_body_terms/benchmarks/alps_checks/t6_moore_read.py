"""Bosonic Moore-Read state from a 3-body (6-operator) NSITETERM.

Bosons in the lowest Landau level on a sphere with 2S flux quanta. The orbitals
m = -S..S are the ALPS sites. The 3-body contact interaction projects every
triple of bosons onto total angular momentum 3S (all three at one point):

    H = sum_M T_M^dag T_M,
    T_M = sum_{m1+m2+m3=M} <S m1; S m2; S m3 | 3S M> b_m1 b_m2 b_m3,
    <S m1; S m2; S m3 | 3S M> = sqrt(C(2S,S+m1) C(2S,S+m2) C(2S,S+m3) / C(6S,3S+M)).

The whole Hamiltonian is one NSITETERM on one n-site group holding every
orbital; ALPS expands the products of sums into strings of 6 boson operators,
several of which can act on the same orbital.

Checks:
1. The full ALPS spectrum equals a dense numpy reference that applies the T_M
   directly to Fock states.
2. The number of exact zero-energy states equals the number of (2,2)-admissible
   occupation patterns: no more than 2 bosons in any 2 consecutive orbitals
   (generalized Pauli principle, Bernevig & Haldane, PRL 100, 246802 (2008)).
   At 2S = N - 2 this is 1: the Moore-Read ground state. Larger 2S gives the
   quasihole zero modes.
"""
import itertools
import os
import sys
from math import comb, sqrt

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun

ZERO = 1e-9


def triples(two_s):
    """{M: {sorted orbital triple (site indices): coefficient in T_M}}.

    The coefficient sums the CG coefficient over the orderings of the triple,
    so T_M = sum over sorted triples of coef * b b b."""
    s = two_s / 2
    norb = two_s + 1
    c = [sqrt(comb(two_s, k)) for k in range(norb)]          # site k is m = k - S
    out = {}
    for k1, k2, k3 in itertools.combinations_with_replacement(range(norb), 3):
        big_m = (k1 + k2 + k3) - 3 * s
        mult = len(set(itertools.permutations((k1, k2, k3))))
        coef = mult * c[k1] * c[k2] * c[k3] / sqrt(comb(3 * two_s, round(3 * s + big_m)))
        out.setdefault(big_m, {})[(k1, k2, k3)] = coef
    return out


def nsite_term(two_s):
    terms = []
    for big_m, tr in sorted(triples(two_s).items()):
        ann = " + ".join(f"{w!r}*b(o{a})*b(o{b})*b(o{c})" for (a, b, c), w in tr.items())
        cre = " + ".join(f"{w!r}*bdag(o{a})*bdag(o{b})*bdag(o{c})" for (a, b, c), w in tr.items())
        terms.append(f"({cre})*({ann})")
    return " + ".join(terms)


def alps_spectrum(n, two_s):
    norb = two_s + 1
    sites = "".join(f'<SITE id="{k + 1}"/>' for k in range(norb))
    lattice = (f'<LATTICES><GRAPH name="LLL sphere" vertices="{norb}">'
               f'<NSITE type="0">{sites}</NSITE></GRAPH></LATTICES>')
    roles = " ".join(f"o{k}" for k in range(norb))
    ham = (f'<HAMILTONIAN name="3-body contact"><BASIS ref="boson"/>'
           f'<NSITETERM sites="{roles}">{nsite_term(two_s)}</NSITETERM></HAMILTONIAN>')
    return alpsrun.run(lattice, ham, {"GRAPH": "LLL sphere", "MODEL": "3-body contact",
                                      "Nmax": n, "CONSERVED_QUANTUMNUMBERS": "N", "N_total": n})


def fock_states(n, norb):
    return [occ for occ in itertools.product(range(n + 1), repeat=norb) if sum(occ) == n]


def reference_spectrum(n, two_s):
    norb = two_s + 1
    upper = fock_states(n, norb)
    lower = {occ: i for i, occ in enumerate(fock_states(n - 3, norb))}
    h = np.zeros((len(upper), len(upper)))
    for tr in triples(two_s).values():
        t = np.zeros((len(lower), len(upper)))                # T_M: N -> N-3
        for col, occ in enumerate(upper):
            for orbs, w in tr.items():
                o, amp = list(occ), w
                for k in reversed(orbs):
                    if o[k] == 0:
                        amp = 0
                        break
                    amp *= sqrt(o[k])
                    o[k] -= 1
                if amp:
                    t[lower[tuple(o)], col] += amp
        h += t.T @ t
    return np.linalg.eigvalsh(h)


def admissible(n, norb):
    return sum(all(occ[k] + occ[k + 1] <= 2 for k in range(norb - 1))
               for occ in fock_states(n, norb))


def main():
    ok = True
    print(f"{'N':>2} {'2S':>3} {'dim':>5} {'max|E_alps-E_ref|':>18} {'zero modes':>11} "
          f"{'admissible':>11} {'E_0':>10} {'gap':>10}")
    for n, two_s in [(4, 2), (4, 3), (4, 4), (6, 4), (6, 5), (6, 6), (8, 6)]:
        e = np.array(alps_spectrum(n, two_s))
        ref = reference_spectrum(n, two_s)
        diff = np.max(np.abs(e - ref)) if len(e) == len(ref) else float("inf")
        zeros = int(np.sum(np.abs(e) < ZERO))
        count = admissible(n, two_s + 1)
        gap = e[zeros] if zeros < len(e) else float("nan")
        good = diff < 1e-9 and zeros == count and (two_s != n - 2 or zeros == 1)
        ok &= good
        print(f"{n:>2} {two_s:>3} {len(e):>5} {diff:>18.2e} {zeros:>11} {count:>11} "
              f"{e[0]:>10.2e} {gap:>10.6f}{'' if good else '   <-- FAIL'}")
    print("ALL OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
