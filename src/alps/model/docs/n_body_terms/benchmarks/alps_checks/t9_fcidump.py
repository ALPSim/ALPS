"""alps2fcidump: ALPS models exported as FCIDUMP, solved by PySCF FCI.

    python t9_fcidump.py          # needs pyscf

Each case writes an ALPS model (site, bond and n-site terms) to an FCIDUMP
with alps2fcidump, diagonalises it with PySCF's FCI (pyscf.fci.direct_spin1,
which fixes (N_alpha, N_beta) like ALPS), and compares the lowest levels with
ALPS fulldiag in the same sector. That fulldiag run sets FCIDUMP="file",
which must write the same file. Then models that FCIDUMP cannot hold
(spins, a Zeeman field, a spin-dependent or a three-body n-site term, and
singlet pair hopping, whose (pq|rs) lacks the 8-fold symmetry of real
orbitals) must be refused, by alps2fcidump and by fulldiag with FCIDUMP set.
"""
import itertools
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun                                          # noqa: E402
from pyscf import fci                                   # noqa: E402
from pyscf.tools import fcidump                         # noqa: E402

NROOTS = 4

HOP = "-t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"
PAIR_K = """-K*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))
            + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )"""
# hopping across the n-site group on one rung pair, correlated with hopping on the
# other: (E_il + E_li)(E_jk + E_kj), spin-independent, a genuine 4-site term
E = "(cdag_up({0})*c_up({1}) + cdag_up({1})*c_up({0}) + cdag_down({0})*c_down({1}) + cdag_down({1})*c_down({0}))"
CORR = "K*" + E.format("l", "i") + "*" + E.format("k", "j")
# the walkthrough's term: up electrons on i, l and down on j, k only
DEMO_K = "K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + cdag_up(i)*cdag_down(j)*c_down(k)*c_up(l))"


def ladder_model(kterm, site="U*n_up(i)*n_down(i)"):
    return f"""<HAMILTONIAN name="fd ladder">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="0"/>
  <PARAMETER name="K" default="0"/><PARAMETER name="h" default="0"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">{site}</SITETERM>
  <BONDTERM source="i" target="j">{HOP}</BONDTERM>
  <NSITETERM sites="i j k l">{kterm}</NSITETERM>
</HAMILTONIAN>"""


def random_nsite_group(m, order, seed):
    """Random spin-independent two-body term on one m-site group whose
    roles sit on the sites in `order`, plus random hopping on a chain."""
    rng = np.random.default_rng(seed)
    roles = [f"r{k}" for k in range(m)]
    g = rng.normal(size=(m,) * 4)
    g = g + g.transpose(1, 0, 2, 3); g = g + g.transpose(0, 1, 3, 2); g = g + g.transpose(2, 3, 0, 1)
    out = []
    for p, q, r, s in itertools.product(range(m), repeat=4):
        for a in ("up", "down"):
            for b in ("up", "down"):
                out.append(f"({float(0.5 * g[p, q, r, s])!r})*cdag_{a}({roles[p]})*cdag_{b}({roles[r]})"
                           f"*c_{b}({roles[s]})*c_{a}({roles[q]})")
    sites = "".join(f'<SITE id="{s + 1}"/>' for s in order)
    edges = "".join(f'<EDGE type="0" source="{k + 1}" target="{k + 2}"/>' for k in range(m - 1))
    lat = (f'<LATTICES><GRAPH name="fd graph" vertices="{m}">{edges}'
           f'<NSITE type="0">{sites}</NSITE></GRAPH></LATTICES>')
    ham = (f'<HAMILTONIAN name="fd random"><PARAMETER name="t" default="1"/><BASIS ref="fermion"/>'
           f'<SITETERM site="i">-0.3*n(i)</SITETERM><BONDTERM source="i" target="j">{HOP}</BONDTERM>'
           f'<NSITETERM sites="{" ".join(roles)}">{" + ".join(out)}</NSITETERM></HAMILTONIAN>')
    return lat, ham


def alps2fcidump(lattice_xml, hamiltonian_xml, params):
    """Run alps2fcidump; returns (FCIDUMP text or None, stderr)."""
    tmp = tempfile.mkdtemp(prefix="fd_")
    try:
        p = {}
        if lattice_xml is not None:
            with open(os.path.join(tmp, "lattices.xml"), "w") as f:
                f.write(lattice_xml)
            p["LATTICE_LIBRARY"] = "lattices.xml"
        if hamiltonian_xml is not None:
            with open(os.path.join(tmp, "models.xml"), "w") as f:
                f.write(alpsrun.model_library(hamiltonian_xml))
            p["MODEL_LIBRARY"] = "models.xml"
        p.update(params)
        with open(os.path.join(tmp, "parm"), "w") as f:
            for k, v in p.items():
                f.write(f'{k}="{v}"\n' if isinstance(v, str) else f"{k}={v!r}\n")
        r = subprocess.run([os.path.join(alpsrun.BIN, "alps2fcidump"), "parm", "FCIDUMP"],
                           cwd=tmp, capture_output=True, text=True)
        if r.returncode:
            return None, r.stderr.strip()
        return os.path.join(tmp, "FCIDUMP"), tmp
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


def pyscf_levels(path):
    d = fcidump.read(path, verbose=False)
    n, ms2 = d["NORB"], d["MS2"]
    na, nb = (d["NELEC"] + ms2) // 2, (d["NELEC"] - ms2) // 2
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-12
    e, _ = solver.kernel(d["H1"], d["H2"], n, (na, nb), ecore=d["ECORE"], nroots=NROOTS)
    return np.sort(np.atleast_1d(e))


CHAIN = '<LATTICES><GRAPH name="fd chain" vertices="6">' + "".join(
    f'<EDGE type="0" source="{k}" target="{k + 1}"/>' for k in range(1, 6)) + "</GRAPH></LATTICES>"

cases = [  # label, lattice xml, hamiltonian xml (None: stock models.xml), params
    ("Hubbard chain L=6 (stock), (3,3)", CHAIN, None,
     {"GRAPH": "fd chain", "MODEL": "fermion Hubbard", "t": 1.0, "U": 4.0, "mu": 0.3,
      "Nup_total": 3, "Ndown_total": 3}),
    ("Hubbard chain L=6 (stock), (3,1)", CHAIN, None,
     {"GRAPH": "fd chain", "MODEL": "fermion Hubbard", "t": 1.0, "U": 4.0, "mu": 0.3,
      "Nup_total": 3, "Ndown_total": 1}),
    ("2x3 ladder, correlated hop K=0.7, (2,2)", alpsrun.ladder_lattice(), ladder_model(CORR),
     {"LATTICE": "ms-ladder", "L": 3, "t": 1.0, "U": 8.0, "K": 0.7, "Nup_total": 2, "Ndown_total": 2}),
    ("2x3 periodic, correlated hop K=0.9, (3,2)", alpsrun.ladder_lattice(boundary="periodic"),
     ladder_model(CORR),
     {"LATTICE": "ms-ladder", "L": 3, "t": 1.0, "U": 8.0, "K": 0.9, "Nup_total": 3, "Ndown_total": 2}),
    ("2x4 ladder, correlated hop K=-1.1, (3,3)", alpsrun.ladder_lattice(), ladder_model(CORR),
     {"LATTICE": "ms-ladder", "L": 4, "t": 1.0, "U": 8.0, "K": -1.1, "Nup_total": 3, "Ndown_total": 3}),
    ("random 3-site (pq|rs), roles 2,0,1, (2,1)", *random_nsite_group(3, [2, 0, 1], 11),
     {"GRAPH": "fd graph", "Nup_total": 2, "Ndown_total": 1}),
    ("random 4-site (pq|rs), roles 3,1,0,2, (2,2)", *random_nsite_group(4, [3, 1, 0, 2], 12),
     {"GRAPH": "fd graph", "Nup_total": 2, "Ndown_total": 2}),
]

refused = [  # label, lattice, hamiltonian, params, expected text in the error
    ("Heisenberg chain (spin basis)", CHAIN, None,
     {"GRAPH": "fd chain", "MODEL": "spin", "J": 1.0}, "fermion basis"),
    ("Zeeman field h*(n_up - n_down)", alpsrun.ladder_lattice(),
     ladder_model(PAIR_K, "U*n_up(i)*n_down(i) - h*(n_up(i) - n_down(i))"),
     {"LATTICE": "ms-ladder", "L": 2, "U": 4.0, "h": 0.3}, "differs between up and down"),
    # SU(2) symmetric but (pq|rs) != (qp|rs): FCIDUMP's 8-fold symmetry fails
    ("singlet pair hopping (t1_ladder)", alpsrun.ladder_lattice(), ladder_model(PAIR_K),
     {"LATTICE": "ms-ladder", "L": 3, "U": 8.0, "K": 0.7}, "symmetry (pq|rs) = (qp|rs)"),
    ("walkthrough term (spin-dependent)", alpsrun.ladder_lattice(), ladder_model(DEMO_K),
     {"LATTICE": "ms-ladder", "L": 2, "U": 4.0, "K": 0.5}, "spin-independent"),
    ("three-body n(i)*n(j)*n(k)", alpsrun.ladder_lattice(), ladder_model("K*n(i)*n(j)*n(k)"),
     {"LATTICE": "ms-ladder", "L": 2, "U": 4.0, "K": 0.5}, "two-body"),
]

bad = 0
for label, lat, ham, params in cases:
    p = dict(params)
    if ham is not None:
        p["MODEL"] = ham.split('name="')[1].split('"')[0]
    path, tmp = alps2fcidump(lat, ham, p)
    if path is None:
        print(f"{label:46} alps2fcidump failed: {tmp}  <-- FAIL")
        bad += 1
        continue
    ref = pyscf_levels(path)
    # the same model through fulldiag's FCIDUMP parameter
    param_file = os.path.join(tmp, "FCIDUMP.param")
    run = dict(p, CONSERVED_QUANTUMNUMBERS="Nup,Ndown", FCIDUMP=param_file)
    run.pop("LATTICE_LIBRARY", None)
    a = np.asarray(alpsrun.run(lat, ham, run)[:len(ref)])
    same = os.path.exists(param_file) and open(param_file).read() == open(path).read()
    shutil.rmtree(tmp, ignore_errors=True)
    d = np.max(np.abs(a - ref))
    bad += d > 1e-9 or not same
    print(f"{label:46} E0 ALPS {a[0]:15.10f}  PySCF FCI {ref[0]:15.10f}  "
          f"max|dE| ({len(ref)} levels) {d:.1e}" + ("" if d < 1e-9 else "  <-- FAIL")
          + ("  FCIDUMP= same" if same else "  FCIDUMP= differs  <-- FAIL"))

for label, lat, ham, params, expect in refused:
    p = dict(params)
    if ham is not None:
        p["MODEL"] = ham.split('name="')[1].split('"')[0]
    path, msg = alps2fcidump(lat, ham, p)
    ok = path is None and expect in msg
    if path is not None:
        shutil.rmtree(msg, ignore_errors=True)
    bad += not ok
    print(f"refused: {label:37} " + ("yes" if ok else "NO  <-- FAIL") + f"   [{(msg or '')[:110]}]")
    # fulldiag with FCIDUMP set must stop with the same reason, writing nothing
    tmp = tempfile.mkdtemp(prefix="fd_")
    param_file = os.path.join(tmp, "FCIDUMP")
    run = dict(p, FCIDUMP=param_file)
    run.pop("LATTICE_LIBRARY", None)
    try:
        alpsrun.run(lat, ham, run)
        err = None
    except RuntimeError as e:
        err = str(e)
    ok = err is not None and "FCIDUMP: " in err and expect in err and not os.path.exists(param_file)
    shutil.rmtree(tmp, ignore_errors=True)
    bad += not ok
    print(f"   and fulldiag FCIDUMP= {'':23}" + ("yes" if ok else "NO  <-- FAIL"))

print("ALL OK" if not bad else f"{bad} FAILED")
sys.exit(1 if bad else 0)
