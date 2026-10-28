# n-body (n-site) terms: benchmarks

Supporting material for `<NSITE>` / `<NSITETERM>` (syntax:
[`../nsite_terms.md`](../nsite_terms.md)). Nothing here is built or
installed; the feature itself is the code under `src/alps/lattice` and
`src/alps/model`. It is opt-in: inputs without `<NSITE>` and
`<NSITETERM>` behave exactly as before.

- `benchmarks/`
  - `casci/compare_casci_nsite.py`: molecules (H3+, H4, LiH, H6), full
    (pq|rs), ALPS `fulldiag` vs QDK CASCI vs a dense Jordan-Wigner matrix.
    `compare_casci_dmrg.py`: the same molecules with ALPS `dmrg`.
    `compare_casci_mps.py`: the CASCI wave function vs the MPS that `dmrg`
    writes with `SAVE_MPS=1`, coefficient by coefficient (needs a `dmrg`
    built with ALPSim/ALPS#216).
    `qdk2alps.py` writes the ALPS input from a QDK Hamiltonian, and
    `qdk_alps_ed.py` registers ALPS as a QDK solver. Needs `qdk-chemistry`.
  - `alps_checks/t1_ladder.py`: ladders with a pair-hopping n-site term vs
    the independent prototype `hyper_ed.py`.
  - `alps_checks/t2_analytic.py`: 4-site term with a known exact spectrum, and
    a particle-number-changing term in the full Fock space.
  - `alps_checks/t5_interleaved.py`: random two-body terms on 2- and 3-site
    n-site groups vs a dense Jordan-Wigner matrix.
  - `alps_checks/t6_moore_read.py`: bosonic Moore-Read state from a 3-body
    (6-operator) contact interaction in the lowest Landau level on a sphere.
    Full spectrum vs a dense reference; zero-mode counts vs the
    (2,2)-admissible patterns, with the unique Moore-Read ground state at
    2S = N - 2. N = 4, 6, 8; takes about 15 seconds.
  - `alps_checks/t7_pairing.py`: pairing terms as site, bond and n-site group
    terms on an open graph and a periodic ladder (with and without momentum
    sectors), full spectrum vs a dense Jordan-Wigner matrix; with Nup, Ndown
    conserved every one of them is refused.
  - `alps_checks/t9_fcidump.py`: `alps2fcidump` (ALPS model -> FCIDUMP, see
    [`../nsite_terms.md`](../nsite_terms.md)). Hubbard chains,
    ladders with a 4-site correlated-hopping n-site term and random
    3- and 4-site (pq|rs) n-site groups: ALPS `fulldiag` vs PySCF FCI on the
    exported FCIDUMP; `fulldiag` with `FCIDUMP="file"` must write the same
    file. Models FCIDUMP cannot hold must be refused by both. Needs
    `pyscf`; takes about a minute.
  - `results/fcidump_solvers.txt`: exported FCIDUMP files solved by PySCF,
    block2, QuSpin and XDiag, compared with ALPS `fulldiag`/`sparsediag`.
    A second table: the same files, and Hubbard rings L=6 and L=10 at U=4,
    mapped to qubit Hamiltonians by Jordan-Wigner in Qiskit Nature and
    OpenFermion, agree with ALPS to 3e-13; on the rings, ALPS's ground-state
    energy matches the Bethe ansatz (Lieb-Wu equations) to 1e-13.
  - `alps_checks/t7_dmrg.py`: n-site terms in `dmrg`: ladders vs
    `fulldiag`, random 4- and 6-operator fermion terms on scrambled n-site groups
    vs a dense matrix, Moore-Read, a 2x6 ladder vs `sparsediag`, and the
    refusal and cutoff warning.
  - `alps_checks/t8_mps.py`: ground-state wave functions from `dmrg`'s MPS
    (`SAVE_MPS=1`, needs ALPSim/ALPS#216) vs exact states built in ALPS's
    local basis (`mpsref.py`): fermions (Hubbard ladders with a pair-hopping
    term, random 3-body terms), bosons (Moore-Read, Bose-Hubbard ring hop) and
    spins (S=1/2 ring exchange, S=1 three-site term).
  - `results/`: output of the runs above.

Running: build ALPS with applications (`cmake --preset default
-DALPS_BUILD_APPLICATIONS=ON`). The scripts use `_build/default/bin` by
default; set `ALPS_BIN` to point elsewhere. `ALPS_XML_PATH` must point to
`src/alps/resources`.
