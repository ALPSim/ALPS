"""Random Hermitian two-body terms on 2- and 3-site plaquettes, all index patterns."""
import itertools, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun
from refjw import dense

rng = np.random.default_rng(7)
bad = 0
for m, order in ((2, [0, 1]), (2, [1, 0]), (3, [0, 1, 2]), (3, [2, 0, 1])):
    roles = [f"r{k}" for k in range(m)]
    g = rng.normal(size=(m,) * 4)
    g = g + g.transpose(1, 0, 2, 3); g = g + g.transpose(0, 1, 3, 2); g = g + g.transpose(2, 3, 0, 1)
    out = []
    for p, q, r, s in itertools.product(range(m), repeat=4):
        for a in ("up", "down"):
            for b in ("up", "down"):
                out.append(f"({float(0.5 * g[p, q, r, s])!r})*cdag_{a}({roles[p]})*cdag_{b}({roles[r]})*c_{b}({roles[s]})*c_{a}({roles[q]})")
    expr = " + ".join(out)
    lat = (f'<LATTICES><GRAPH name="gm" vertices="{m}">'
           f'<PLAQUETTE type="0">{"".join(f"<SITE id={chr(34)}{s + 1}{chr(34)}/>" for s in order)}</PLAQUETTE></GRAPH></LATTICES>')
    ham = f'<HAMILTONIAN name="x"><BASIS ref="fermion"/><PLAQUETTETERM sites="{" ".join(roles)}">{expr}</PLAQUETTETERM></HAMILTONIAN>'
    a = alpsrun.run(lat, ham, {"GRAPH": "gm", "MODEL": "x", "CONSERVED_QUANTUMNUMBERS": ""})
    ref = dense(m, [], [order], "0", expr, roles, {})
    d = np.max(np.abs(np.array(a) - ref))
    bad += d > 1e-9
    print(f"m={m} roles at {order}: levels {len(a)}  max|dE| {d:.1e}")
print("ALL OK" if not bad else f"{bad} FAILED")
