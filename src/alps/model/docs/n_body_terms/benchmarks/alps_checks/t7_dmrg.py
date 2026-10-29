"""NSITETERMs in ALPS dmrg vs exact references.

1. Ladders with a 4-site plaquette term (t1 cases) vs ALPS fulldiag.
2. Random 4-operator and 6-operator (3-body) fermion terms on 4- and 3-site
   n-site groups whose role order is scrambled relative to the DMRG site order,
   plus hopping, vs a dense Jordan-Wigner matrix: ground and first excited
   state in several (Nup, Ndown) sectors, including odd ones.
3. Bosonic Moore-Read (t6): one n-site group over all orbitals, products of 6
   boson operators with several on one orbital. Zero mode and gap vs dense.
4. A 2x6 ladder, where DMRG truncates, vs ALPS sparsediag.
5. dmrg refuses a term that changes N, conserved or not; fulldiag warns.

DMRG needs at least 4 sites (a limit of ALPS dmrg itself, also for plain
SITETERM/BONDTERM models), so 3-site systems are not tested.
"""
import os
import sys
import time

import numpy as np
import scipy.sparse as sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun
import t6_moore_read as mr
from refjw import Op, fock_ops

# tight Lanczos tolerance: the default leaves excited states (E1) good to only ~1e-6
DMRG = {"SWEEPS": 4, "MAXSTATES": 300, "LANCZOS_TOLERANCE": 1e-12}
bad = 0


def check(label, d, tol, extra=""):
    global bad
    ok = d < tol
    bad += not ok
    print(f"  {label:44} max|dE| {d:.1e}{extra}" + ("" if ok else "  <-- FAIL"), flush=True)


# 1. ladders ------------------------------------------------------------------
HOP = "-t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"
DEMO_K = "K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + cdag_up(i)*cdag_down(j)*c_down(k)*c_up(l))"
PAIR_K = """-K*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))
            + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )"""


def ladder_ham(kterm):
    return f"""<HAMILTONIAN name="ms test">
  <PARAMETER name="t" default="1"/><PARAMETER name="U" default="0"/><PARAMETER name="K" default="0"/>
  <BASIS ref="fermion"/>
  <SITETERM site="i">U*n_up(i)*n_down(i)</SITETERM>
  <BONDTERM source="i" target="j">{HOP}</BONDTERM>
  <NSITETERM sites="i j k l">{kterm}</NSITETERM>
</HAMILTONIAN>"""


print("1. ladders, dmrg vs fulldiag (ground state)")
for label, L, bc, U, K, kt, nu, nd in [
        ("2x2 demo K=0.5 (1,1)", 2, "open", 4.0, 0.5, DEMO_K, 1, 1),
        ("2x3 open U=8 K=0.7 (2,2)", 3, "open", 8.0, 0.7, PAIR_K, 2, 2),
        ("2x3 open U=8 K=0.7 (3,2)", 3, "open", 8.0, 0.7, DEMO_K, 3, 2),
        ("2x3 periodic U=8 K=0.9 (3,2)", 3, "periodic", 8.0, 0.9, PAIR_K, 3, 2),
        ("2x4 open U=8 K=1.1 (3,2)", 4, "open", 8.0, 1.1, PAIR_K, 3, 2),
        ("2x4 open U=8 K=1.1 (3,3)", 4, "open", 8.0, 1.1, DEMO_K, 3, 3)]:
    p = {"LATTICE": "ms-ladder", "MODEL": "ms test", "L": L, "U": U, "K": K, "t": 1.0,
         "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", "Nup_total": nu, "Ndown_total": nd}
    lat = alpsrun.ladder_lattice(boundary=bc)
    ed = alpsrun.run(lat, ladder_ham(kt), p)
    dm = alpsrun.run(lat, ladder_ham(kt), dict(p, **DMRG), program="dmrg")
    check(label, abs(dm[0] - ed[0]), 1e-6, f"  E0 {dm[0]:.10f}")

# 2. random 4- and 6-operator terms, scrambled n-site groups ---------------------
print("2. random 2-body and 3-body fermion terms vs dense JW (E0, E1)")
rng = np.random.default_rng(11)
NS = 6


def hermitian(terms):
    out = []
    for w, ops in terms:
        s = "*".join(f"{'cdag' if d else 'c'}_{sp}({r})" for d, sp, r in ops)
        hc = "*".join(f"{'c' if d else 'cdag'}_{sp}({r})" for d, sp, r in reversed(ops))
        out.append(f"({w!r})*({s} + {hc})")
    return " + ".join(out)


R4, R3 = ["a", "b", "c", "d"], ["x", "y", "z"]
t4, t3 = [], []
for _ in range(40):
    p, q, r, s = rng.integers(0, 4, 4)
    a, b = rng.choice(["up", "down"], 2)
    t4.append((float(rng.normal()), [(1, a, R4[p]), (1, b, R4[r]), (0, b, R4[s]), (0, a, R4[q])]))
for _ in range(20):
    k = rng.integers(0, 3, 6)
    sp = rng.choice(["up", "down"], 3)
    t3.append((float(rng.normal()), [(1, sp[0], R3[k[0]]), (1, sp[1], R3[k[1]]), (1, sp[2], R3[k[2]]),
                                     (0, sp[2], R3[k[3]]), (0, sp[1], R3[k[4]]), (0, sp[0], R3[k[5]])]))
E4, E3 = hermitian(t4), hermitian(t3)
HOP1 = "-(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"
P4 = [[3, 0, 5, 1], [2, 4, 1, 5], [5, 4, 3, 2], [0, 1, 2, 3]]
P3 = [[4, 0, 2], [1, 5, 3], [5, 0, 3]]


def site(i):
    return f'<SITE id="{i + 1}"/>'


lat = (f'<LATTICES><GRAPH name="g" vertices="{NS}">'
       + "".join(f'<EDGE type="0" source="{i + 1}" target="{i + 2}"/>' for i in range(NS - 1))
       + "".join(f'<NSITE type="0">{"".join(map(site, p))}</NSITE>' for p in P4)
       + "".join(f'<NSITE type="1">{"".join(map(site, p))}</NSITE>' for p in P3)
       + '</GRAPH></LATTICES>')
ham = (f'<HAMILTONIAN name="x"><BASIS ref="fermion"/><BONDTERM source="i" target="j">{HOP1}</BONDTERM>'
       f'<NSITETERM type="0" sites="{" ".join(R4)}">{E4}</NSITETERM>'
       f'<NSITETERM type="1" sites="{" ".join(R3)}">{E3}</NSITETERM></HAMILTONIAN>')
ops = fock_ops(NS, sparse=True)
H = Op(sparse.csr_matrix((1 << 2 * NS,) * 2))
for i in range(NS - 1):
    H = H + eval(HOP1, {}, {**ops, "i": i, "j": i + 1})
for p in P4:
    H = H + eval(E4, {}, {**ops, **dict(zip(R4, p))})
for p in P3:
    H = H + eval(E3, {}, {**ops, **dict(zip(R3, p))})
assert abs(H.m - H.m.T).max() < 1e-12
UP, DN = int("01" * NS, 2), int("10" * NS, 2)
for nu, nd in [(1, 1), (2, 1), (3, 2), (3, 3), (4, 2)]:
    st = [s for s in range(1 << 2 * NS) if bin(s & UP).count("1") == nu and bin(s & DN).count("1") == nd]
    ref = np.linalg.eigvalsh(H.m[st][:, st].toarray())
    dm = alpsrun.run(lat, ham, {"GRAPH": "g", "MODEL": "x", "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown",
                                "Nup_total": nu, "Ndown_total": nd, "NUMBER_EIGENVALUES": 2, **DMRG},
                     program="dmrg")
    check(f"6 sites, sector ({nu},{nd})", max(abs(dm[0] - ref[0]), abs(dm[1] - ref[1])), 1e-6,
          f"  E0 {dm[0]:.10f}")

# 3. bosonic Moore-Read --------------------------------------------------------
print("3. bosonic Moore-Read, dmrg vs dense (E0, gap)")
for n, two_s in [(6, 4), (8, 6)]:
    norb = two_s + 1
    lat = (f'<LATTICES><GRAPH name="LLL sphere" vertices="{norb}">'
           f'<NSITE type="0">{"".join(map(site, range(norb)))}</NSITE></GRAPH></LATTICES>')
    ham = (f'<HAMILTONIAN name="3-body contact"><BASIS ref="boson"/>'
           f'<NSITETERM sites="{" ".join(f"o{k}" for k in range(norb))}">{mr.nsite_term(two_s)}'
           f'</NSITETERM></HAMILTONIAN>')
    ref = mr.reference_spectrum(n, two_s)
    dm = alpsrun.run(lat, ham, {"GRAPH": "LLL sphere", "MODEL": "3-body contact", "Nmax": n,
                                "CONSERVED_QUANTUMNUMBERS": "N", "N_total": n, "NUMBER_EIGENVALUES": 2,
                                "SWEEPS": 4, "MAXSTATES": 200, "LANCZOS_TOLERANCE": 1e-12}, program="dmrg")
    check(f"N={n} 2S={two_s} (Moore-Read ground state)", max(abs(dm[0] - ref[0]), abs(dm[1] - ref[1])), 1e-5,
          f"  E0 {dm[0]:.1e} gap {dm[1]:.6f}")

# 4. truncation matters --------------------------------------------------------
print("4. 2x6 ladder, dmrg vs sparsediag")
p = {"LATTICE": "ms-ladder", "MODEL": "ms test", "L": 6, "U": 8.0, "K": 1.1, "t": 1.0,
     "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", "Nup_total": 4, "Ndown_total": 4}
lat = alpsrun.ladder_lattice()
t0 = time.time()
ed = alpsrun.run(lat, ladder_ham(PAIR_K), p, program="sparsediag")
t_ed = time.time() - t0
t0 = time.time()
dm = alpsrun.run(lat, ladder_ham(PAIR_K), dict(p, SWEEPS=6, MAXSTATES=400), program="dmrg")
check("2x6 open U=8 K=1.1 (4,4), 400 states", abs(dm[0] - ed[0]), 1e-5,
      f"  E0 {dm[0]:.10f}  sparsediag {t_ed:.0f}s dmrg {time.time() - t0:.0f}s")

# 5. refusals and warnings -----------------------------------------------------
print("5. refusals and warnings")
lat4 = f'<LATTICES><GRAPH name="g4" vertices="4"><NSITE type="0">{"".join(map(site, range(4)))}</NSITE></GRAPH></LATTICES>'
pair = '<HAMILTONIAN name="b"><BASIS ref="boson"/><NSITETERM sites="i j k l">bdag(i)*bdag(j)*b(k) + bdag(k)*b(j)*b(i)</NSITETERM></HAMILTONIAN>'
for label, conserved in (("conserved", {"CONSERVED_QUANTUMNUMBERS": "N", "N_total": 2}),
                         ("not conserved", {"CONSERVED_QUANTUMNUMBERS": ""})):
    try:
        alpsrun.run(lat4, pair, {"GRAPH": "g4", "MODEL": "b", "Nmax": 2, **conserved, **DMRG}, program="dmrg")
        msg = ""
    except RuntimeError as e:
        msg = str(e)
    ok = "changes the total N" in msg
    bad += not ok
    print(f"  {'N-changing term, N ' + label + ', refused':44} {'yes' if ok else 'NO  <-- FAIL'}")


def stderr_of(program):
    import subprocess
    real = subprocess.run

    def capture(*a, **k):
        r = real(*a, **k)
        captured.append(r.stderr)
        return r
    captured = []
    subprocess.run = capture
    try:
        alpsrun.run(lat4, pair, {"GRAPH": "g4", "MODEL": "b", "Nmax": 2, "CONSERVED_QUANTUMNUMBERS": "",
                                 "NUMBER_EIGENVALUES": 1, **DMRG}, program=program)
    finally:
        subprocess.run = real
    return "".join(captured)


ok = "Warning: an <NSITETERM> changes the total N" in stderr_of("fulldiag")
bad += not ok
print(f"  {'Nmax cutoff warning in fulldiag':44} {'yes' if ok else 'NO  <-- FAIL'}")

print("ALL OK" if not bad else f"{bad} FAILED")
sys.exit(1 if bad else 0)
