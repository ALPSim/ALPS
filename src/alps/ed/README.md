# `ALPS::ed` — exact diagonalization of compiled operators

Exact diagonalization for models built with [`ALPS::models`](../models/README.md).
The component takes a `CompiledOperator`/`CompiledModel` (numeric local
factors with Jordan–Wigner strings), enumerates a product basis, a charge
sector or a charge-and-momentum-projected basis, assembles the sparse matrix
and solves it densely (LAPACK) or iteratively (Lanczos over the in-tree
IETL). Spectra can be resolved by every conserved charge and lattice momentum
at once, and thermodynamic averages follow from the levels. It does not
depend on the XML catalog, legacy parameters, HDF5 or the simulation runtime.
The Python package binds it as `pyalps.ed`, with NumPy and SciPy views of
the results ([pyalps README](../../../python/pyalps/README.md#models-exact-diagonalization-and-classic-dmrg)).

```c++
#include <alps/ed.hpp>
#include <alps/models.hpp>

namespace lat = alps::lattice; namespace mdl = alps::models;
auto const model = mdl::heisenberg(lat::chain(16));          // periodic, S=1/2
alps::ed::ProductOptions sector; sector.sector = {{"Sz", 0.}};
auto const problem = alps::ed::prepare_product(model, sector); // dim 12870
alps::ed::EigenstateOptions request; request.count = 3;
auto const result = problem.eigenstates(request);           // Lanczos above dim 1000
double const e0 = result.eigenvalues[0];
auto const sz_sq = alps::ed::expectation_values(result, model.operator_decomposition(my_observable).op);

// The same sector at momentum pi (about 800 states), with eigenvectors in
// the projected basis, and the full charge/momentum-resolved spectrum.
alps::ed::ProjectedOptions at_pi; at_pi.sector = {{"Sz", 0.}}; at_pi.momentum = {std::acos(-1.)};
auto const block = alps::ed::prepare_projected(model, at_pi).eigenstates(request);
alps::ed::ModelSpectrumOptions all; all.all_momenta = true; all.count = 4;
auto const spectrum = alps::ed::model_spectrum(model, all);   // blocks labelled by Sz and k
auto const thermal = alps::ed::thermal_surface(levels, charges, temperatures, fields);
```

## Headers (`<alps/ed/...>`, namespace `alps::ed`)

| header | content |
| --- | --- |
| `product_basis.hpp` | `enumerate_product_basis(site_dims, constraints)`: full product basis or a charge sector, enumerated with pruning; site 0 is the most significant digit |
| `sparse_from_local.hpp` | `visit_local_transitions`, `visit_local_operator`, `build_coo_from_local` (any operator on an explicit basis, complex COO) and `build_sparse_from_local` (real CSR from a `CompiledModel` with QN constraints) |
| `sparse_matrix.hpp` | `SparseMatrix<Scalar>`, `canonical_csr` (sorted, summed, zero-free CSR), `build_csr_from_local`, `build_canonical_csr` (the same CSR assembled row by row from the operator's adjoint, without the uncombined COO list; product problems use it) |
| `eigensolver.hpp` | `EigenstateOptions`, `dense_eigensystem` (LAPACK `dsyevd`/`zheevd`), `select_eigenvalues` |
| `lanczos.hpp` | `lanczos_csr`: the deflated IETL Lanczos driver (real symmetric and complex Hermitian CSR) |
| `product_problem.hpp` | `prepare_product` (`Hamiltonian` or `TermSum` + sites), `ProductProblem::eigenstates`, `operator_matrix`, `expectation_values` |
| `spectrum.hpp` | `SpectrumLocalSpace`, `SpectrumTranslation`, `SpectrumBlock` (momentum, charges, method and, for point-group blocks, the irrep label, dimension and group order); `symmetry_spectrum`: charge- and momentum-resolved eigenvalues of a product-basis CSR matrix |
| `finite_group.hpp` | `FiniteGroupIrrep`: the closure of a faithful site-permutation group with one unitary irrep, checked for generator relations, unitarity and irreducibility; `dihedral_irrep` (D_n on a ring), `validate_irrep_family` (complete, inequivalent family of one group) |
| `projected_basis.hpp` | `ProjectedBasis`/`make_projected_basis` (orbit representatives of a charge sector under commuting translations at one momentum), `make_finite_group_basis` (the full isotypic component of one irrep), `projected_operator`, `projected_sparse_operator`, `apply_projected_operator`, `expect_projected_operator`, `expand_projected_state`, `symmetry_spectrum_from_local`, `finite_spectrum_from_local` |
| `projected_problem.hpp` | `select_momentum`, `projected_spaces`, `prepare_projected` and `prepare_finite` (`Hamiltonian`, `TermSum` + sites, or a compiled operator on `ProjectedSpaces` or an explicit basis), `ProjectedProblem::eigenstates`/`hamiltonian`, `model_spectrum`, `finite_spectrum` |
| `thermodynamics.hpp` | header-only `thermal_weights` and `thermal_surface` (energy, free energy, entropy, specific heat, charge and susceptibility on a temperature × field grid, signed temperatures), `ThermalAccumulator` |
| `ed.hpp` | umbrella |

## Conventions

- Matrices are stored with row = output state and column = input state;
  `canonical_csr` sorts columns within a row, sums duplicates in their
  contribution (term) order and drops exact zeros; `build_canonical_csr`
  gives bit-identical offsets, indices and values.
- `enumerate_product_basis`, `ProductProblem` and `ProjectedBasis` use site 0
  as the most significant digit; `build_sparse_from_local` keeps the kit's
  least-significant-first mixed-radix codes (its states are returned).
- Fermionic terms carry their Jordan–Wigner string sites; the parity
  diagonal of each site's basis supplies the signs, exactly as in the
  kit's dense reference assembly (`src/alps/models/tests/dense_e0.hpp`).
  Translations of graded local spaces pick up the permutation sign of the
  occupied odd states.
- Eigenvalues are ascending; eigenvectors are column-major
  (`vectors[column * dimension + row]`). Projected eigenvectors are complex
  for every momentum.
- Translations map source site to translated site and must commute; a
  momentum is a row of the Cartesian grid `k_d = 2π j_d / L_d` (last
  direction fastest), the translation eigenvalue is `exp(+i k)`, and
  `select_momentum` snaps a wavevector in radians onto that grid. A
  `Hamiltonian` supplies the generators of its periodic lattice directions.
- Operators on a projected basis are validated before assembly: Hermiticity,
  conservation of every declared charge (forbidden amplitudes are summed
  exactly in rational arithmetic, so cancellations between huge terms are
  neither trusted nor fabricated) and invariance under the translations.
- A finite point group is given by faithful site permutations (source site
  to image site) and one unitary irrep whose row-major matrices satisfy
  `D(gh) = D(g) D(h)`; `FiniteGroupIrrep` closes the group and checks the
  relations, unitarity and the character norm. `make_finite_group_basis`
  keeps every row of the isotypic component (an orbit-local central
  projector, its range resolved by twice-reorthogonalized Gram–Schmidt), so
  multiplets stay explicit and `finite_spectrum` over a complete family
  (`validate_irrep_family`) reproduces the whole sector. Invariance under the
  generators is checked exactly in rational arithmetic: a scalar offset or a
  cancellation between huge coefficients cannot hide a symmetry-breaking
  term.
- `eigenstates` chooses Lanczos automatically above dimension 1000 and
  only when the dimension exceeds `8*count + 200`; `which` is `SA`, `LA`,
  `LM` or `SM` (dense only). Requested Lanczos solves on smaller spaces are
  dense and say so in the result. Spectrum blocks (`symmetry_spectrum`,
  `symmetry_spectrum_from_local`, `model_spectrum`) follow the same automatic
  policy and report the method per block.

## The Lanczos driver

`lanczos_csr` wraps `ietl::lanczos` (the sparsediag engine) in a deflation
loop. One IETL run converges the lowest distinct eigenvalues and their Ritz
vectors; converged vectors are locked, the next run works on the projected,
shifted operator `(1-P)(±A-σ)(1-P)`, and a final run on the remaining space
verifies that no lower eigenvalue was missed. This resolves degenerate
multiplets, which a single Krylov sequence cannot see, at the cost of a few
extra runs. Start vectors are deterministic, so repeated solves agree
bit for bit. `LA` runs on `-A`; `LM` combines both ends; `SM` is refused
(no shift-invert).

## Origin and provenance

Ported from the private skilledwolf/alpsx fork (branch
`agent/v3-execution-and-api-audit`, commit f8ee6e98b, directory
`src/solvers/diag/`); every ported file names its source. Adaptations: C++17
(no `std::span`, defaulted comparisons, templated lambdas or designated
initializers), `<alps/ed/...>` include names, `ALPS_ED_DECL` exports; the
dense solver uses LAPACK only (alpsx used Eigen for small blocks and
`zheevr`), the Lanczos driver is new (alpsx used Spectra over Eigen,
which this tree does not vendor), projected blocks are held as canonical
CSR instead of Eigen sparse matrices, and the finite-group closure and
projector ranges use small local dense helpers (alpsx used Eigen). The
catalog importer (`<alps/model/catalog.hpp>`, in the legacy model module)
supplies compiled operators, parities and charges of XML catalog models; see
"Catalog models" below. Not ported: the alpsx diag runtime's catalog problem
adapters and HDF5 persistence.

## Catalog models

`alps::models::CatalogModel<>` (header `<alps/model/catalog.hpp>`, library
`ALPS::alps`) resolves `MODEL`/`LATTICE` parameters through the XML catalogs
and compiles the Hamiltonian, or any catalog operator, into a
`CompiledOperator`. Its `local_parity()`, `conserved_quantum_numbers()`,
`local_charge(name)` and `sector()` fill a `ProjectedSpaces` and a sector for
`prepare_projected`, `model_spectrum`, `prepare_finite` and the CSR entry
points, so a catalog model and a kit `Hamiltonian` share every solver:

```cpp
alps::models::CatalogModel<> model(parameters);   // MODEL, LATTICE, L, Sz_total, ...
auto const op = model.compile_hamiltonian();
alps::ed::ProjectedSpaces spaces;
spaces.charges = model.conserved_quantum_numbers();
auto const parity = model.local_parity();
for (std::size_t site = 0; site < parity.size(); ++site) {
    alps::ed::SpectrumLocalSpace space;
    space.parity = parity[site];
    for (auto const& name : spaces.charges) space.quantum_numbers.push_back(model.local_charge(name)[site]);
    spaces.local_spaces.push_back(space);
}
alps::ed::ProjectedOptions options;
options.sector = model.sector();
auto const problem = alps::ed::prepare_projected(op, spaces, {}, options);
```

The same compiled operators and physical bases drive the classic DMRG
solver ([`ALPS::classic_dmrg`](../classic_dmrg/README.md)) for chains and
ladders beyond the reach of exact diagonalization. A catalog model can also
run on a kit lattice: `CatalogModel<>(lattice, parameters)` takes an
`alps::lattice::Lattice` in place of the catalog's `LATTICE`
(`<alps/lattice/kit.hpp>` describes the bridge).

The legacy module's tests lock the importer against the kit factories
(`catalog_model_parity`), against `hamiltonian_matrix` element by element on
the same catalog basis (`ed_vs_hamiltonian_matrix`), and the kit lattices
against the catalog graphs (`lattice_catalog_structure`).

## Tests

`ctest -L ed`: ported unit tests (basis enumeration, local transitions,
streamed operators, canonical CSR, term accumulation, thermodynamics), the
sparse-assembly lock against the dense Kronecker reference, dense and Lanczos
solver locks (analytic chains, degenerate multiplets, selection and
validation), product problems with observables, the end-to-end lock
reproducing the ground state energies of the 14 frozen model fixtures, and
the projected-basis suite: fermionic rings against analytic dispersions,
reflection and two-direction translation groups, conserved-sector escape
detection (exact cancellation versus hidden residues), momentum-resolved
spectra of the CSR and local-operator paths against each other and against
the unprojected sector, lattice-generator agreement for honeycomb and ladder
geometries, and the dispatch policy of the IETL driver; and the finite-group
suite: dihedral, cyclic and tetrahedral irreps against independent central
projectors, complete families reproducing the unprojected spectra,
cross-irrep observables, exact rejection of broken symmetries and a
polarized 100-site sector.
