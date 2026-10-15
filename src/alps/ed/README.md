# `ALPS::ed` — exact diagonalization of compiled operators

Exact diagonalization for models built with [`ALPS::models`](../models/README.md).
The component takes a `CompiledOperator`/`CompiledModel` (numeric local
factors with Jordan–Wigner strings), enumerates a product basis or a charge
sector, assembles the sparse matrix and solves it densely (LAPACK) or
iteratively (Lanczos over the in-tree IETL). It does not depend on the XML
catalog, legacy parameters, HDF5 or the simulation runtime.

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
```

## Headers (`<alps/ed/...>`, namespace `alps::ed`)

| header | content |
| --- | --- |
| `product_basis.hpp` | `enumerate_product_basis(site_dims, constraints)`: full product basis or a charge sector, enumerated with pruning; site 0 is the most significant digit |
| `sparse_from_local.hpp` | `visit_local_transitions`, `visit_local_operator`, `build_coo_from_local` (any operator on an explicit basis, complex COO) and `build_sparse_from_local` (real CSR from a `CompiledModel` with QN constraints) |
| `sparse_matrix.hpp` | `SparseMatrix<Scalar>`, `canonical_csr` (sorted, summed, zero-free CSR), `build_csr_from_local` |
| `eigensolver.hpp` | `EigenstateOptions`, `dense_eigensystem` (LAPACK `dsyevd`/`zheevd`), `select_eigenvalues` |
| `lanczos.hpp` | `lanczos_csr`: the deflated IETL Lanczos driver (real symmetric and complex Hermitian CSR) |
| `product_problem.hpp` | `prepare_product` (`Hamiltonian` or `TermSum` + sites), `ProductProblem::eigenstates`, `operator_matrix`, `expectation_values` |
| `ed.hpp` | umbrella |

## Conventions

- Matrices are stored with row = output state and column = input state;
  `canonical_csr` sorts columns within a row, sums duplicates and drops
  exact zeros.
- `enumerate_product_basis` and `ProductProblem` use site 0 as the most
  significant digit; `build_sparse_from_local` keeps the kit's
  least-significant-first mixed-radix codes (its states are returned).
- Fermionic terms carry their Jordan–Wigner string sites; the parity
  diagonal of each site's basis supplies the signs, exactly as in the
  kit's dense reference assembly (`src/alps/models/tests/dense_e0.hpp`).
- Eigenvalues are ascending; eigenvectors are column-major
  (`vectors[column * dimension + row]`).
- `eigenstates` chooses Lanczos automatically above dimension 1000 and
  only when the dimension exceeds `8*count + 200`; `which` is `SA`, `LA`,
  `LM` or `SM` (dense only). Requested Lanczos solves on smaller spaces are
  dense and say so in `ProductEigenResult::method`.

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
`zheevr`) and the Lanczos driver is new (alpsx used Spectra over Eigen,
which this tree does not vendor). Not ported yet: momentum/point-group
projected bases and spectra, thermodynamics, the catalog (XML) problem
adapters and HDF5 persistence.

## Tests

`ctest -L ed`: ported unit tests (basis enumeration, local transitions,
streamed operators, canonical CSR, term accumulation), the sparse-assembly
lock against the dense Kronecker reference, dense and Lanczos solver locks
(analytic chains, degenerate multiplets, selection and validation), product
problems with observables, and the end-to-end lock reproducing the ground
state energies of the 14 frozen model fixtures.
