"""Ground-state wave functions of NSITETERM models: dmrg's MPS vs exact.

dmrg writes its final ground state as an MPS with SAVE_MPS=1 (ALPSim/ALPS#216),
so this needs a dmrg built with that change as well as the n-site terms;
set ALPS_BIN to it. Each case contracts the MPS to a full state vector and
compares it with the exact ground state, built independently in the same
local basis (mpsref.py): |<exact|MPS>| must be 1 and every component must
agree, which checks the sign of every amplitude, not only the energy.

Fermions, bosons and spins, each with an n-site term:
- Hubbard ladders with a 4-site pair-hopping term (2x3, 2x4).
- Random 4-operator and 6-operator fermion terms on scrambled n-site groups.
- Bosonic Moore-Read: 6-operator 3-body term on all orbitals.
- Bose-Hubbard ladder with a 4-site ring hop, Nmax = 2.
- Spin-1/2 ladder with 4-spin ring exchange.
- Spin-1 chain with a 3-site (S_i.S_j)(S_j.S_k) term.

Molecules against the QDK CASCI wave function: ../casci/compare_casci_mps.py.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import alpsrun
import mpsref
import t6_moore_read as mr

DMRG = {"SWEEPS": 4, "MAXSTATES": 300, "LANCZOS_TOLERANCE": 1e-12}
TOL = 1e-8
bad = 0
print(f"{'case':46}{'dim':>7}{'|<exact|MPS>|':>16}{'max|dpsi|':>11}{'<H>_MPS-E0':>12}{'E_dmrg-E0':>11}")


def check(label, lattice, ham, params, nsites, terms, sector, S=None):
    """Run dmrg, rebuild H in the MPS's local basis, compare ground states."""
    global bad
    energies, mps = alpsrun.run_mps(lattice, ham, dict(DMRG, **params))
    types = set(mps["site_type"])
    assert len(types) == 1, "one site type expected"
    basis = mps["site_basis"][types.pop()]
    H = mpsref.hamiltonian(nsites, basis, terms, S=S)
    e, psi = mpsref.ground_state(H, nsites, basis, sector)
    assert e[1] - e[0] > 1e-6, f"{label}: degenerate ground state"
    ov, emps, dev, _ = mpsref.compare(mps, psi, H)
    # components converge as the square root of the DMRG energy error; a wrong
    # sign anywhere moves them by ~0.1
    ok = abs(1 - ov) < TOL and dev < 1e-4 and abs(emps - e[0]) < TOL and abs(energies[0] - e[0]) < TOL
    bad += not ok
    print(f"{label:46}{len(e):>7}{ov:16.12f}{dev:11.1e}"
          f"{emps - e[0]:12.1e}{energies[0] - e[0]:11.1e}" + ("" if ok else "  <-- FAIL"), flush=True)


def site(i):
    return f'<SITE id="{i + 1}"/>'


def graph(name, n, edges, groups):
    """GRAPH with edges [(i, j, type)] and n-site groups {type: [site lists]}."""
    return (f'<LATTICES><GRAPH name="{name}" vertices="{n}">'
            + "".join(f'<EDGE type="{t}" source="{i + 1}" target="{j + 1}"/>' for i, j, t in edges)
            + "".join(f'<NSITE type="{t}">{"".join(map(site, p))}</NSITE>'
                      for t, ps in groups.items() for p in ps)
            + '</GRAPH></LATTICES>')


def ladder(L):
    """2xL ladder, site 2x+y; legs, rungs and the L-1 square plaquettes (i j k l around)."""
    edges = [(2 * x, 2 * x + 1, 1) for x in range(L)]
    edges += [(2 * x + y, 2 * x + 2 + y, 0) for x in range(L - 1) for y in (0, 1)]
    plaq = [[2 * x, 2 * x + 1, 2 * x + 3, 2 * x + 2] for x in range(L - 1)]
    return edges, plaq


# --- fermions: Hubbard ladders with pair hopping ------------------------------
HOP = "-t*(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"
PAIR_K = ("-K*( (cdag_up(l)*cdag_down(k) - cdag_down(l)*cdag_up(k)) * (c_down(j)*c_up(i) - c_up(j)*c_down(i))"
          " + (cdag_up(i)*cdag_down(j) - cdag_down(i)*cdag_up(j)) * (c_down(k)*c_up(l) - c_up(k)*c_down(l)) )")
DEMO_K = "K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + cdag_up(i)*cdag_down(j)*c_down(k)*c_up(l))"
for L, kterm, kname, U, K, nu, nd in [(3, PAIR_K, "pair", 8.0, 0.7, 2, 2), (3, DEMO_K, "demo", 4.0, 0.5, 2, 1),
                                      (4, PAIR_K, "pair", 8.0, 1.1, 3, 3)]:
    edges, plaq = ladder(L)
    ham = (f'<HAMILTONIAN name="f"><PARAMETER name="t" default="1"/><BASIS ref="fermion"/>'
           f'<SITETERM site="i">U*n_up(i)*n_down(i)</SITETERM>'
           f'<BONDTERM source="i" target="j">{HOP}</BONDTERM>'
           f'<NSITETERM sites="i j k l">{kterm}</NSITETERM></HAMILTONIAN>')
    p = {"t": 1.0, "U": U, "K": K}
    terms = [(f"{U!r}*n_up(i)*n_down(i)", "i", [(s,) for s in range(2 * L)]),
             (HOP.replace("t*", "1.0*"), "ij", [(i, j) for i, j, _ in edges]),
             (kterm.replace("K*", f"{K!r}*"), "ijkl", plaq)]
    check(f"fermions: 2x{L} Hubbard + {kname} term ({nu},{nd})", graph("lad", 2 * L, edges, {0: plaq}), ham,
          dict(p, GRAPH="lad", MODEL="f", CONSERVED_QUANTUMNUMBERS="Nup,Ndown", Nup_total=nu, Ndown_total=nd),
          2 * L, terms, {"Nup": nu, "Ndown": nd})

# --- fermions: random 2-body and 3-body terms, scrambled n-site groups -----------
rng = np.random.default_rng(5)


def hermitian(terms):
    out = []
    for w, ops in terms:
        s = "*".join(f"{'cdag' if d else 'c'}_{sp}({r})" for d, sp, r in ops)
        hc = "*".join(f"{'c' if d else 'cdag'}_{sp}({r})" for d, sp, r in reversed(ops))
        out.append(f"({w!r})*({s} + {hc})")
    return " + ".join(out)


R4, R3 = "abcd", "xyz"
t4, t3 = [], []
for _ in range(30):
    a, b, c, d = rng.integers(0, 4, 4)
    s1, s2 = rng.choice(["up", "down"], 2)
    t4.append((float(rng.normal()), [(1, s1, R4[a]), (1, s2, R4[c]), (0, s2, R4[d]), (0, s1, R4[b])]))
for _ in range(15):
    k = rng.integers(0, 3, 6)
    s = rng.choice(["up", "down"], 3)
    t3.append((float(rng.normal()), [(1, s[0], R3[k[0]]), (1, s[1], R3[k[1]]), (1, s[2], R3[k[2]]),
                                     (0, s[2], R3[k[3]]), (0, s[1], R3[k[4]]), (0, s[0], R3[k[5]])]))
E4, E3 = hermitian(t4), hermitian(t3)
HOP1 = "-(cdag_up(i)*c_up(j) + cdag_up(j)*c_up(i) + cdag_down(i)*c_down(j) + cdag_down(j)*c_down(i))"
P4, P3 = [[3, 0, 5, 1], [2, 4, 1, 5], [0, 1, 2, 3]], [[4, 0, 2], [5, 1, 3]]
chain = [(i, i + 1, 0) for i in range(5)]
ham = (f'<HAMILTONIAN name="r"><BASIS ref="fermion"/><BONDTERM source="i" target="j">{HOP1}</BONDTERM>'
       f'<NSITETERM type="0" sites="a b c d">{E4}</NSITETERM>'
       f'<NSITETERM type="1" sites="x y z">{E3}</NSITETERM></HAMILTONIAN>')
terms = [(HOP1, "ij", [(i, j) for i, j, _ in chain]), (E4, R4, P4), (E3, R3, P3)]
for nu, nd in [(2, 1), (3, 2), (3, 3)]:
    check(f"fermions: random 4-op + 6-op terms ({nu},{nd})", graph("r", 6, chain, {0: P4, 1: P3}), ham,
          {"GRAPH": "r", "MODEL": "r", "CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", "Nup_total": nu, "Ndown_total": nd},
          6, terms, {"Nup": nu, "Ndown": nd})

# --- bosons: Moore-Read --------------------------------------------------------
for n, two_s in [(6, 4), (8, 6)]:
    norb = two_s + 1
    roles = [f"o{k}" for k in range(norb)]
    ham = (f'<HAMILTONIAN name="mr"><BASIS ref="boson"/>'
           f'<NSITETERM sites="{" ".join(roles)}">{mr.nsite_term(two_s)}</NSITETERM></HAMILTONIAN>')
    check(f"bosons: Moore-Read N={n} 2S={two_s} (3-body, 6 ops)", graph("s", norb, [], {0: [list(range(norb))]}),
          ham, {"GRAPH": "s", "MODEL": "mr", "Nmax": n, "CONSERVED_QUANTUMNUMBERS": "N", "N_total": n},
          norb, [(mr.nsite_term(two_s), roles, [tuple(range(norb))])], {"N": n})

# --- bosons: Bose-Hubbard ladder with a ring hop --------------------------------
RING_B = "K*(bdag(i)*b(j)*bdag(k)*b(l) + bdag(j)*b(i)*bdag(l)*b(k))"
BHOP = "-t*(bdag(i)*b(j) + bdag(j)*b(i))"
edges, plaq = ladder(3)
ham = (f'<HAMILTONIAN name="bh"><BASIS ref="boson"/><SITETERM site="i">0.5*U*n(i)*(n(i)-1)</SITETERM>'
       f'<BONDTERM source="i" target="j">{BHOP}</BONDTERM>'
       f'<NSITETERM sites="i j k l">{RING_B}</NSITETERM></HAMILTONIAN>')
terms = [("0.5*3.0*n(i)*(n(i)-1)", "i", [(s,) for s in range(6)]),
         (BHOP.replace("t*", "1.0*"), "ij", [(i, j) for i, j, _ in edges]),
         (RING_B.replace("K*", "0.8*"), "ijkl", plaq)]
for n in (3, 5):
    check(f"bosons: 2x3 Bose-Hubbard + ring hop, N={n}, Nmax=2", graph("bl", 6, edges, {0: plaq}), ham,
          {"GRAPH": "bl", "MODEL": "bh", "t": 1.0, "U": 3.0, "K": 0.8, "Nmax": 2,
           "CONSERVED_QUANTUMNUMBERS": "N", "N_total": n}, 6, terms, {"N": n})

# --- spins ---------------------------------------------------------------------
def dot(a, b):
    return f"(Sz({a})*Sz({b}) + 0.5*(Splus({a})*Sminus({b}) + Sminus({a})*Splus({b})))"


RING_S = f"K*({dot('i', 'j')}*{dot('k', 'l')} + {dot('i', 'l')}*{dot('j', 'k')} - {dot('i', 'k')}*{dot('j', 'l')})"
edges, plaq = ladder(4)
ham = (f'<HAMILTONIAN name="sr"><BASIS ref="spin"/>'
       f'<BONDTERM source="i" target="j">J*{dot("i", "j")}</BONDTERM>'
       f'<NSITETERM sites="i j k l">{RING_S}</NSITETERM></HAMILTONIAN>')
terms = [(dot("i", "j"), "ij", [(i, j) for i, j, _ in edges]), (RING_S.replace("K*", "0.6*"), "ijkl", plaq)]
check("spins: S=1/2 2x4 ladder + 4-spin ring exchange", graph("sl", 8, edges, {0: plaq}), ham,
      {"GRAPH": "sl", "MODEL": "sr", "local_S": 0.5, "J": 1.0, "K": 0.6,
       "CONSERVED_QUANTUMNUMBERS": "Sz", "Sz_total": 0}, 8, terms, {"Sz": 0})

THREE = f"K*({dot('i', 'j')}*{dot('j', 'k')} + {dot('j', 'k')}*{dot('i', 'j')})"
chain = [(i, i + 1, 0) for i in range(5)]
triples = [[i, i + 1, i + 2] for i in range(4)]
ham = (f'<HAMILTONIAN name="s1"><BASIS ref="spin"/>'
       f'<BONDTERM source="i" target="j">J*{dot("i", "j")}</BONDTERM>'
       f'<NSITETERM sites="i j k">{THREE}</NSITETERM></HAMILTONIAN>')
terms = [(dot("i", "j"), "ij", [(i, j) for i, j, _ in chain]), (THREE.replace("K*", "0.3*"), "ijk", triples)]
check("spins: S=1 chain L=6 + 3-site (Si.Sj)(Sj.Sk)", graph("c1", 6, chain, {0: triples}), ham,
      {"GRAPH": "c1", "MODEL": "s1", "local_S": 1, "J": 1.0, "K": 0.3,
       "CONSERVED_QUANTUMNUMBERS": "Sz", "Sz_total": 0}, 6, terms, {"Sz": 0}, S=1)

print("ALL OK" if not bad else f"{bad} FAILED")
sys.exit(1 if bad else 0)
