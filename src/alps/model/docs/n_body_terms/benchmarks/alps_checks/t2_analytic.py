"""Analytic and independent-reference checks of ALPS PLAQUETTETERM.

1. Analytic: open 4-site chain, H = -t sum hop + U n_{a,up} n_{a,dn}, where
   a = sum_i phi_i c_i is the k-th hopping eigenmode. Written in the site basis
   the U term is a single 4-site term (ALPS expands the product of sums), with
   pieces on 1, 2, 3 and 4 distinct sites. Exact spectrum:
       E = sum_{occupied up modes} eps + sum_{occupied dn modes} eps + U [a up and dn occupied]
   Run with the plaquette's roles in site order and in a scrambled order.
2. Non-conserving 4-site term (changes N by 4) with CONSERVED_QUANTUMNUMBERS=""
   against a dense Jordan-Wigner matrix of the same strings, and with Nup,Ndown
   conserved, which must be refused.
"""
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun

L = 4
HOP = "-t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"


def graph(order, name="g4"):
    edges = "".join(f'<EDGE type="0" source="{i + 1}" target="{i + 2}"/>' for i in range(L - 1))
    plaq = "".join(f'<SITE id="{s + 1}"/>' for s in order)
    return (f'<LATTICES><GRAPH name="{name}" vertices="{L}">{edges}'
            f'<PLAQUETTE type="0">{plaq}</PLAQUETTE></GRAPH></LATTICES>')


def ham(plaq_term, roles="i j k l"):
    return f"""<HAMILTONIAN name="ms analytic">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="0"/>
  <BASIS ref="fermion"/>
  <BONDTERM source="i" target="j">{HOP}</BONDTERM>
  <PLAQUETTETERM sites="{roles}">{plaq_term}</PLAQUETTETERM>
</HAMILTONIAN>"""


def mode_term(order, phi):
    """U * Adag_up A_up Adag_dn A_dn with A = sum_r phi[site of role r] c(role r)."""
    roles = "ijkl"
    def lin(op, spin):
        return "(" + " + ".join(f"({float(phi[s])!r})*{op}_{spin}({roles[r]})" for r, s in enumerate(order)) + ")"
    return f"U*{lin('cdag', 'up')}*{lin('c', 'up')}*{lin('cdag', 'down')}*{lin('c', 'down')}"


def analytic(eps, ka, U, nu, nd):
    out = []
    for up in itertools.combinations(range(L), nu):
        for dn in itertools.combinations(range(L), nd):
            out.append(sum(eps[k] for k in up) + sum(eps[k] for k in dn) + U * (ka in up and ka in dn))
    return np.sort(out)


# ------------------------------------------------------------- dense reference

from refjw import Op, fock_ops, dense  # noqa: E402


# ------------------------------------------------------------------- run

bad = 0
def report(label, a, ref, tol=1e-9):
    global bad
    d = np.max(np.abs(np.asarray(a) - np.asarray(ref))) if len(a) == len(ref) else np.inf
    bad += not d < tol
    print(f"{label:52} levels {len(a):4d}  E0 {a[0]:14.10f}  max|dE| {d:.1e}" + ("" if d < tol else "  <-- FAIL"))

t, U = 1.0, 3.7
eps = [-2 * t * np.cos(np.pi * (k + 1) / (L + 1)) for k in range(L)]
for ka in (0, 1):
    phi = [np.sqrt(2 / (L + 1)) * np.sin(np.pi * (ka + 1) * (i + 1) / (L + 1)) for i in range(L)]
    for order in ([0, 1, 2, 3], [2, 0, 3, 1], [3, 2, 1, 0]):
        for nu, nd in ((1, 1), (2, 1), (2, 2), (3, 2), (3, 1)):
            a = alpsrun.run(graph(order), ham(mode_term(order, phi)),
                            {"GRAPH": "g4", "MODEL": "ms analytic", "t": t, "U": U,
                             "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", "Nup_total": nu, "Ndown_total": nd})
            report(f"analytic: U on mode k={ka}, roles at {order}, ({nu},{nd})", a, analytic(eps, ka, U, nu, nd))

# 2. non-conserving term, full Fock space
NC = "D*(cdag_up(i)*cdag_down(j)*cdag_up(k)*cdag_down(l) + c_down(l)*c_up(k)*c_down(j)*c_up(i)) + W*(cdag_up(i)*c_up(k)*cdag_down(l)*c_down(j) + cdag_down(j)*c_down(l)*cdag_up(k)*c_up(i))"
for order in ([0, 1, 2, 3], [1, 3, 0, 2]):
    p = {"t": 1.0, "D": 0.8, "W": -0.45}
    a = alpsrun.run(graph(order), ham(NC).replace('<PARAMETER name="U" default="0"/>',
                                                   '<PARAMETER name="D" default="0"/><PARAMETER name="W" default="0"/>'),
                    {"GRAPH": "g4", "MODEL": "ms analytic", **p, "CONSERVED_QUANTUMNUMBERS": ""})
    ref = dense(L, [(i, i + 1) for i in range(L - 1)], [order], HOP, NC, "ijkl", p)
    report(f"non-conserving, no sector, roles at {order}", a, ref)

try:
    alpsrun.run(graph([0, 1, 2, 3]), ham(NC).replace('<PARAMETER name="U" default="0"/>',
                                                       '<PARAMETER name="D" default="0"/><PARAMETER name="W" default="0"/>'),
                {"GRAPH": "g4", "MODEL": "ms analytic", "t": 1.0, "D": 0.8, "W": 0.0,
                 "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", "Nup_total": 1, "Ndown_total": 1})
    print("non-conserving with Nup,Ndown conserved: accepted  <-- FAIL"); bad += 1
except RuntimeError as e:
    msg = [l for l in str(e).splitlines() if "PLAQUETTETERM" in l]
    ok = bool(msg)
    bad += not ok
    print("non-conserving with Nup,Ndown conserved: refused" + (f": {msg[0].strip()[:110]}..." if ok else f"  <-- wrong error: {e}"))

print("ALL OK" if not bad else f"{bad} FAILED")
