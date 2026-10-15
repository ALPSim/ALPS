# `ALPS::classic_dmrg` — classic DMRG on compiled models

Finite-system, two-site density-matrix renormalization group for models built
with [`ALPS::models`](../models/README.md): Adrian Feiguin's dmtk engine (the
one behind the `dmrg` application) driven by the numeric local terms of a
`CompiledOperator` on its `PhysicalBasis`. The engine warms up, sweeps with a
kept-state schedule, targets one or several lowest states of one conserved
sector and measures compiled operators on them. The component depends on
`ALPS::models` and, privately, on Eigen; it does not use the XML catalog,
legacy parameters, HDF5, LAPACK or the simulation runtime, and consumers of
the installed SDK need neither Eigen nor LAPACK.

```c++
#include <alps/dmrg/classic.hpp>
#include <alps/models.hpp>

namespace lat = alps::lattice; namespace mdl = alps::models; namespace cd = alps::dmrg::classic;
auto const model = mdl::hubbard(lat::chain(32, 1.0, lat::Boundary::open), 1.0, 4.0);
cd::Options options;
options.sweeps = cd::ramped_sweeps(6, 400);        // 2 half sweeps per sweep, states ramp to 400
options.sector = {{"N", 32.}, {"Sz", 0.}};
options.num_eigenvalues = 2;                       // ground state and gap
cd::Measurements observe;
observe.terms["n"] = cd::site_terms(model, "n");               // one column per site
observe.terms["SzSz"] = cd::pair_terms(model, "Sz", "Sz");     // one column per ordered pair
observe.expectation["H"] = mdl::compile_hamiltonian(model).op; // <H> per target
auto const result = cd::solve(model, options, observe, [](cd::Iteration const& it) { /* progress */ });
double const e0 = result.energies[0], gap = result.energies[1] - e0;
auto const& density = result.measurements.at("n").values[0];  // row = target state
```

## Interface (`<alps/dmrg/classic.hpp>`, namespace `alps::dmrg::classic`)

| name | content |
| --- | --- |
| `Options` | `sweeps` (`fixed_sweeps`, `ramped_sweeps` or an explicit `SweepStates` list; the right-to-left half sweep runs first), `warmup_states`, `num_eigenvalues`, optional `truncation_error` with `truncation_max_states`, optional `lanczos_tolerance`, `arithmetic` (real or complex), `sector` (name → total of a declared conserved charge), `scratch_parent`, `log`, `verbose` |
| `Measurements` | `expectation`: name → operator, one value `<ψ\|O\|ψ>`; `terms`: name → operator, one value per term (one column per site, pair, …) |
| `site_terms(model, op)`, `pair_terms(model, a, b)` | builders over a model's local operator names: one term per site, one term per ordered site pair (i-major, on-site product at i = j); no two columns are ever combined |
| `solve(op, basis, options, measurements, callback)` | the engine on a compiled operator and its physical basis |
| `solve(model, …)` | validates Hermiticity, compiles the `Hamiltonian` under its declared charges and convention, solves on `physical_basis(model)` |
| `Result` | ascending `energies`, `measurements` (rows follow the targets, columns the terms with their factor `sites`), the `iterations` observed (direction, index, energy, truncation error, entropy), `final_truncation_error`, the resolved `options` |
| `IterationCallback` | observes every completed iteration; may throw to abort the solve |

## What the engine accepts

- At least four sites. Every term acts on at most two distinct sites. A local
  factor is split into components of definite charge change (for the
  sector's charges and the fermionic charge) and definite parity change, so
  `Sx Sx` is fine; the components that leave the sector are dropped from
  measurements and must cancel in the Hamiltonian. Like terms are summed, the
  Hamiltonian must be Hermitian and must conserve every charge of the sector.
- Fermions under the strict convention need an integer charge whose parity is
  the fermion parity on every site (the particle number `N` of the built-in
  sites); it becomes the engine's fermionic quantum number, conserved or not.
  Two odd factors must carry the strict convention's Jordan–Wigner string
  `[i, j)` that `ALPS::models` emits; odd terms cannot be represented. Under
  the lattice convention every operator is bosonic.
- At most four quantum numbers are tracked: the sector's charges plus the
  fermionic charge. Charges and targets are half-integers in `short` range;
  an empty sector or more targets than the sector holds are rejected before
  the engine starts.
- Real arithmetic needs real coefficients and factor matrices (`Sp Sm` rather
  than `Sy Sy`); complex ones select `Arithmetic::complex`.
- Measured values are real parts. A two-site product inside one block is
  measured for a second target from a stored block operator, which the
  solver keeps only when several targets and two-site measurements meet.
- Solves are serialized process-wide (the engine keeps global state) and may
  not be nested; each solve owns a private scratch directory under
  `scratch_parent` (default: the system temporary directory) that is removed
  when it returns, and restores the engine's state afterwards, including
  after an exception.

## Origin and provenance

The engine is `src/alps/classic_dmrg/src/dmtk/`, the Eigen-based rework of
`src/apps/dmrg/dmrg/dmtk` ported from the private skilledwolf/alpsx fork
(branch `agent/v3-execution-and-api-audit`, commit f8ee6e98b, directory
`src/solvers/tensor_network/src/classic_dmrg/`); every ported file names its
source. The driver (`src/classic.cpp`) follows alpsx's `DMRG<T>::run`,
`create_site_operator`, `build_*_operator` and `save_results` call for call,
fed from a `CompiledOperator` instead of `graph_helper`/`model_helper` and
`alps::SiteOperator`; the invocation guard is alpsx's without its expression
stream. Adaptations: `dmtk::half_integer` and `dmtk/support.hpp` replace the
ALPS model and utility headers (no legacy headers are included), the
scratch-directory check uses `std::filesystem`, measurements are compiled
operators (alpsx took catalog expressions), operator matrices are addressed
through the engine's charge-sorted block basis (the legacy site basis was
already sorted, the kit orders spin states from +S down), the stored-product
path measures after the final sweep, and the catalog snapshot, random
disorder state, HDF5 output and run recipe are not ported. The XML `dmrg`
application is unchanged and keeps its own copy of dmtk.

## Tests

`ctest -L classic_dmrg`: energies against `ALPS::ed` with enough kept states
for exact sweeps (open Heisenberg chain and ladder, spin-1 chain with
anisotropy, sector targets, an unrestricted solve, three targeted states,
hard-core bosons, the Hubbard chain at and away from half filling, frustrated
spinless fermions on a next-nearest-neighbour chain, complex hopping phases,
truncated ramped sweeps converging from above, iteration observation and
abort); measurements against `ALPS::ed` (site and pair spin terms, fermion
densities and hopping correlations, split charge components, scalar terms,
one row per target); and the rejected inputs. The SDK consumer contract
solves the open 6-site chain through the installed component.
