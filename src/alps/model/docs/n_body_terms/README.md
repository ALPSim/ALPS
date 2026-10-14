# n-body (plaquette) terms: walkthrough and benchmarks

Supporting material for `<PLAQUETTE>` / `<PLAQUETTETERM>` (syntax:
[`../plaquette_terms.md`](../plaquette_terms.md)). Nothing here is built or
installed; the feature itself is the code under `src/alps/lattice` and
`src/alps/model`. It is opt-in: inputs without `<PLAQUETTE>` and
`<PLAQUETTETERM>` behave exactly as before.

- `ALPS_n_body_terms.pdf`: one 2x2 plaquette with a 4-site quartic term,
  followed through ALPS step by step (input files, local table, matrix fill,
  eigenvalues), then the benchmarks and the list of changes.
- `pdf_source/`: the example input (`lattice.xml`, `model.xml`, `parm`), the
  program that prints ALPS's basis, table and matrix (`dump_plaquette.cpp`),
  its output (`data/`), and `make_pdf.py`, which builds the PDF from them.
- `benchmarks/`
  - `casci/compare_casci_nsite.py`: molecules (H3+, H4, LiH, H6), full
    (pq|rs), ALPS `fulldiag` vs QDK CASCI vs a dense Jordan-Wigner matrix.
    `qdk2alps.py` writes the ALPS input from a QDK Hamiltonian, and
    `qdk_alps_ed.py` registers ALPS as a QDK solver. Needs `qdk-chemistry`.
  - `alps_checks/t1_ladder.py`: ladders with a pair-hopping plaquette term vs
    the independent prototype `hyper_ed.py`.
  - `alps_checks/t2_analytic.py`: 4-site term with a known exact spectrum, and
    a particle-number-changing term in the full Fock space.
  - `alps_checks/t5_interleaved.py`: random two-body terms on 2- and 3-site
    plaquettes vs a dense Jordan-Wigner matrix.
  - `alps_checks/t6_moore_read.py`: bosonic Moore-Read state from a 3-body
    (6-operator) contact interaction in the lowest Landau level on a sphere.
    Full spectrum vs a dense reference; zero-mode counts vs the
    (2,2)-admissible patterns, with the unique Moore-Read ground state at
    2S = N - 2. Takes about 7 minutes.
  - `results/`: output of the runs above.

Running: build ALPS with applications (`cmake --preset default
-DALPS_BUILD_APPLICATIONS=ON`). The scripts use `_build/default/bin` by
default; set `ALPS_BIN` to point elsewhere. `ALPS_XML_PATH` must point to
`src/alps/resources`.
