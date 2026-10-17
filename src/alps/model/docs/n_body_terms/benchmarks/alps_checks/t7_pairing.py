"""Pairing (particle-number-changing) terms, and terms that leave a fixed sector.

1. No fixed sector: pairing as an NSITETERM (2-site singlet, reversed role
   order, on-site pair, a 4-site N+2 product, interleaved roles) and as a
   BONDTERM, on an open 4-site graph and a periodic ladder (L=3), with and
   without momentum sectors. Full spectrum vs a dense Jordan-Wigner matrix.
2. Nup,Ndown conserved: every non-conserving term (site, bond, n-site group; with
   and without momentum sectors) must be refused, not silently projected.
3. Control: a conserving n-site term with Nup,Ndown and momentum sectors;
   the union over all sectors must equal the full dense spectrum.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun
from refjw import Op, fock_ops

HOP = "-t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"
DAG = {"cdag_up": "c_up", "c_up": "cdag_up", "cdag_down": "c_down", "c_down": "cdag_down"}


def with_hc(products):
    """ALPS term string for sum_p coef*prod + h.c."""
    out = []
    for c, prod in products:
        for p in (prod, [(DAG[o], r) for o, r in reversed(prod)]):
            out.append(f"({c!r})*" + "*".join(f"{o}({r})" for o, r in p))
    return " + ".join(out)


NSITE_PAIRING = with_hc([
    (0.7,   [("cdag_up", "i"), ("cdag_down", "j")]),
    (-0.7,  [("cdag_down", "i"), ("cdag_up", "j")]),
    (-0.45, [("cdag_down", "k"), ("cdag_up", "i")]),                                   # reversed role order
    (0.3,   [("cdag_up", "l"), ("cdag_down", "l")]),                                   # on-site pair
    (0.25,  [("cdag_up", "i"), ("cdag_down", "k"), ("cdag_up", "l"), ("c_up", "j")]),  # N+2 on 4 sites
    (0.15,  [("cdag_up", "i"), ("c_up", "j"), ("cdag_down", "i"), ("cdag_down", "k")]),  # interleaved roles
])
NSITE_4SITE = with_hc([(0.25, [("cdag_up", "i"), ("cdag_down", "k"), ("cdag_up", "l"), ("c_up", "j")])])
BOND_PAIRING = with_hc([(0.7, [("cdag_up", "i"), ("cdag_down", "j")]),
                        (-0.7, [("cdag_down", "i"), ("cdag_up", "j")])])
SITE_PAIRING = with_hc([(0.3, [("cdag_up", "i"), ("cdag_down", "i")])])
PAIR_HOP = ("-0.9*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))"
            " + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )")

T, U = 1.0, 2.0


def ham(bond="", group="", site=""):
    return (f'<HAMILTONIAN name="pairing"><PARAMETER name="t" default="1"/><PARAMETER name="U" default="0"/>'
            f'<BASIS ref="fermion"/>'
            f'<SITETERM site="i">U*n_up(i)*n_down(i){" + " + site if site else ""}</SITETERM>'
            f'<BONDTERM source="i" target="j">{HOP}{" + " + bond if bond else ""}</BONDTERM>'
            + (f'<NSITETERM sites="i j k l">{group}</NSITETERM>' if group else "")
            + "</HAMILTONIAN>")


# geometries: (lattice xml, ALPS parameters, nsites, bonds, n-site groups)
OPEN4 = ('<LATTICES><GRAPH name="g4" vertices="4">'
         + "".join(f'<EDGE source="{i + 1}" target="{i + 2}"/>' for i in range(3))
         + '<NSITE><SITE id="1"/><SITE id="2"/><SITE id="3"/><SITE id="4"/></NSITE></GRAPH></LATTICES>',
         {"GRAPH": "g4"}, 4, [(0, 1), (1, 2), (2, 3)], [(0, 1, 2, 3)])
L = 3
s = lambda v, x: 2 * (x % L) + v
LADDER = (alpsrun.ladder_lattice(boundary="periodic"), {"LATTICE": "ms-ladder", "L": L}, 2 * L,
          [(s(0, x), s(1, x)) for x in range(L)] + [(s(v, x), s(v, x + 1)) for v in (0, 1) for x in range(L)],
          [(s(0, x), s(1, x), s(1, x + 1), s(0, x + 1)) for x in range(L)])


ONSITE_U = "U*cdag_up(i)*c_up(i)*cdag_down(i)*c_down(i)"


def dense_full(geo, bond="", group="", site=""):
    """refjw's dense Jordan-Wigner matrix, with site terms added."""
    _, _, n, bonds, groups = geo
    params = {"t": T, "U": U}
    ops = fock_ops(n)
    H = Op(np.zeros((1 << 2 * n,) * 2))
    site_expr = ONSITE_U + (" + " + site if site else "")
    bond_expr = HOP + (" + " + bond if bond else "")
    for i in range(n):
        H = H + eval(site_expr, {}, {**ops, **params, "i": i})
    for i, j in bonds:
        H = H + eval(bond_expr, {}, {**ops, **params, "i": i, "j": j})
    if group:
        for sites in groups:
            H = H + eval(group, {}, {**ops, **params, **dict(zip("ijkl", sites))})
    assert np.allclose(H.m, H.m.T)
    return np.linalg.eigvalsh(H.m)


bad = 0


def report(label, a, ref, tol=1e-9):
    global bad
    d = np.max(np.abs(np.asarray(a) - np.asarray(ref))) if len(a) == len(ref) else np.inf
    bad += not d < tol
    print(f"{label:74} levels {len(a):5d}  E0 {a[0]:14.10f}  max|dE| {d:.1e}" + ("" if d < tol else "  <-- FAIL"))


def refused(label, geo, extra, want, **terms):
    global bad
    try:
        alpsrun.run(geo[0], ham(**terms), {**geo[1], "MODEL": "pairing", "t": T, "U": U,
                                           "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", **extra})
        print(f"{label:74} accepted  <-- FAIL"); bad += 1
    except RuntimeError as e:
        msg = [l.strip() for l in str(e).splitlines() if want in l]
        bad += not msg
        print(f"{label:74} refused: " + (msg[0][:100] + "..." if msg else f"<-- wrong error: {str(e)[-200:]}"))


NOTRANS = {"TRANSLATION_SYMMETRY": "false"}

# 1. no fixed sector
for name, geo, extra in (("open 4-site graph", OPEN4, {}),
                         ("periodic ladder, momentum sectors", LADDER, {}),
                         ("periodic ladder, no translation", LADDER, NOTRANS)):
    for kind, terms in (("NSITETERM pairing", {"group": NSITE_PAIRING}),
                        ("BONDTERM pairing", {"bond": BOND_PAIRING}),
                        ("SITETERM pairing", {"site": SITE_PAIRING})):
        a = alpsrun.run(geo[0], ham(**terms), {**geo[1], "MODEL": "pairing", "t": T, "U": U,
                                               "CONSERVED_QUANTUMNUMBERS": "", **extra})
        report(f"{kind}, {name}", a, dense_full(geo, **terms))

# 2. Nup,Ndown conserved: refused
for name, geo, extra in (("open 4-site graph", OPEN4, {}),
                         ("periodic ladder, momentum sectors", LADDER, {}),
                         ("periodic ladder, no translation", LADDER, NOTRANS)):
    refused(f"NSITETERM pairing, Nup,Ndown, {name}", geo, extra, "<NSITETERM>", group=NSITE_PAIRING)
    refused(f"NSITETERM 4-site N+2 only, Nup,Ndown, {name}", geo, extra, "<NSITETERM>", group=NSITE_4SITE)
    refused(f"BONDTERM pairing, Nup,Ndown, {name}", geo, extra, "<BONDTERM>", bond=BOND_PAIRING)
    refused(f"SITETERM pairing, Nup,Ndown, {name}", geo, extra, "<SITETERM>", site=SITE_PAIRING)

# 3. conserving control: all sectors and momenta together give the full spectrum
a = alpsrun.run(LADDER[0], ham(group=PAIR_HOP), {**LADDER[1], "MODEL": "pairing", "t": T, "U": U,
                                                "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown"})
report("conserving pair hopping, all Nup,Ndown and momentum sectors", a, dense_full(LADDER, group=PAIR_HOP))

print("ALL OK" if not bad else f"{bad} FAILED")
