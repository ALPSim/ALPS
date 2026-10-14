"""ALPS (new PLAQUETTETERM) vs the independent Python prototype, ladder cases."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun
import hyper_ed as ed

HOP = "-t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"
DEMO_K = "K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + cdag_up(i)*cdag_down(j)*c_down(k)*c_up(l))"
PAIR_K = """-K*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))
            + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )"""

def ham(kterm):
    return f"""<HAMILTONIAN name="ms test">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="0"/><PARAMETER name="K" default="0"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">U*n_up(i)*n_down(i)</SITETERM>
  <BONDTERM source="i" target="j">{HOP}</BONDTERM>
  <PLAQUETTETERM sites="i j k l">{kterm}</PLAQUETTETERM>
</HAMILTONIAN>"""

def proto(L, boundary, U, K, kterm, nup, ndn):
    m = {"lattice": {"kind": "ladder", "L": L, "boundary": boundary},
         "sector": {"Nup": nup, "Ndown": ndn},
         "parameters": {"t": 1.0, "U": U, "K": K},
         "term": [{"on": "site", "sites": ["i"], "expr": "U*n_up(i)*n_down(i)"},
                  {"on": "bond", "sites": ["i", "j"], "expr": HOP},
                  {"on": "plaquette", "sites": ["i", "j", "k", "l"], "expr": kterm}]}
    return ed.solve(m, n_levels=200, verbose=False)["energies"]

cases = [  # label, L, boundary, U, K, kterm, nup, ndn, reference E0 from the plan
    ("2x2 plaquette demo K=0.5", 2, "open", 4.0, 0.5, DEMO_K, 1, 1, -3.3597011665),
    ("2x2 plaquette demo K=0",   2, "open", 4.0, 0.0, DEMO_K, 1, 1, -3.4185507189),
    ("2x3 open, U=8, K=0.7",     3, "open", 8.0, 0.7, PAIR_K, 2, 2, -5.3669606044),
    ("2x3 open, U=8, K=-1.3",    3, "open", 8.0, -1.3, PAIR_K, 2, 2, -5.0797988157),
    ("2x3 periodic, U=8, K=0.9", 3, "periodic", 8.0, 0.9, PAIR_K, 3, 2, -5.0117188128),
    ("2x4 open, U=8, K=0",       4, "open", 8.0, 0.0, PAIR_K, 3, 3, -5.7500660284),
    # (d) 3-electron sectors: spectator electron sets signs
    ("2x2 demo K=0.5 (2,1)",     2, "open", 4.0, 0.5, DEMO_K, 2, 1, None),
    ("2x3 open K=0.7 (2,1)",     3, "open", 8.0, 0.7, PAIR_K, 2, 1, None),
    ("2x3 open K=0.7 (3,2)",     3, "open", 8.0, 0.7, DEMO_K, 3, 2, None),
    ("2x3 periodic K=0.9 (2,2)", 3, "periodic", 8.0, 0.9, DEMO_K, 2, 2, None),
    ("2x4 open K=1.1 (3,2)",     4, "open", 8.0, 1.1, PAIR_K, 3, 2, None),
]

bad = 0
for label, L, bc, U, K, kt, nu, nd, ref in cases:
    a = alpsrun.run(alpsrun.ladder_lattice(boundary=bc), ham(kt),
                    {"LATTICE": "ms-ladder", "MODEL": "ms test", "L": L, "U": U, "K": K, "t": 1.0,
                     "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", "Nup_total": nu, "Ndown_total": nd})
    p = sorted(proto(L, bc, U, K, kt, nu, nd))
    n = min(len(a), len(p))
    dspec = max(abs(x - y) for x, y in zip(a[:n], p[:n]))
    dref = abs(a[0] - ref) if ref is not None else 0.0
    ok = dspec < 1e-9 and dref < 1e-9 and len(a) >= len(p) or (len(p) == 200 and dspec < 1e-9)
    bad += not ok
    print(f"{label:28} E0 ALPS {a[0]:15.10f}  proto {p[0]:15.10f}  levels {n:4d}  max|dE| {dspec:.1e}"
          + (f"  vs plan {dref:.1e}" if ref is not None else "") + ("" if ok else "  <-- FAIL"))
print("ALL OK" if not bad else f"{bad} FAILED")
