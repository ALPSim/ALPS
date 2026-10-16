# `ALPS::classic_dmrg` — classic DMRG on compiled models

Finite-system, two-site density-matrix renormalization group for models built
with [`ALPS::models`](../models/README.md): Adrian Feiguin's dmtk engine
driven by the numeric local terms of a `CompiledOperator` on its
`PhysicalBasis`. The XML `dmrg` application runs on this component
([below](#the-dmrg-application)). The engine warms up, sweeps with a
kept-state schedule, targets one or several lowest states of one conserved
sector and measures compiled operators on them. The component depends on
`ALPS::models` and, privately, on Eigen (dense eigensolvers) and the platform
BLAS (block products); it does not use the XML catalog, legacy parameters,
HDF5, LAPACK or the simulation runtime, and consumers of the installed SDK
need neither Eigen nor BLAS headers.
The Python package binds it as `pyalps.classic_dmrg` when the SDK has the
component ([pyalps README](../../../python/pyalps/README.md#models-exact-diagonalization-and-classic-dmrg)).

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
| `Options` | `sweeps` (`fixed_sweeps`, `ramped_sweeps` or an explicit `SweepStates` list; the right-to-left half sweep runs first), `warmup_states`, `num_eigenvalues`, optional `truncation_error` with `truncation_max_states`, optional `lanczos_tolerance` ([below](#the-lanczos-stop-rule)), `arithmetic` (real or complex), `sector` (name → total of a declared conserved charge), `scratch_parent`, `log`, `verbose` |
| `Measurements` | `expectation`: name → operator, one value `<ψ\|O\|ψ>`; `terms`: name → operator, one value per term (one column per site, pair, …) |
| `site_terms(model, op)`, `pair_terms(model, a, b)` | builders over a model's local operator names: one term per site, one term per ordered site pair (i-major, on-site product at i = j); no two columns are ever combined |
| `solve(op, basis, options, measurements, callback)` | the engine on a compiled operator and its physical basis |
| `solve(model, …)` | validates Hermiticity, compiles the `Hamiltonian` under its declared charges and convention, solves on `physical_basis(model)` |
| `Result` | ascending `energies`, `measurements` (rows follow the targets, columns the terms with their factor `sites`), the `iterations` observed (direction, index, energy, truncation error, entropy), `final_truncation_error`, the resolved `options` |
| `IterationCallback` | observes every completed iteration; may throw to abort the solve |

## What the engine accepts

- At least four sites. Every term acts on at most two distinct sites, its
  factors in either site order. A local factor is split into components of
  definite charge change (for the sector's charges and the fermionic
  charges) and definite parity change, so `Sx Sx` is fine; the components
  that leave the sector are dropped from measurements and must cancel in the
  Hamiltonian.
- Like terms are combined before the engine sees them, also terms that
  differ only by a scale moved between coefficient and factor matrices and
  a product written as one term per matrix element of its first factor (the
  XML catalog's bond terms). The sum must be Hermitian, which is checked on
  the dense matrix of every site and site pair, and must conserve every
  charge of the sector.
- Fermions under the strict convention need integer charges whose summed
  parity is the fermion parity on every site: the particle number `N` of
  the built-in sites, or else the smallest such set, such as the XML
  catalog's `Nup` and `Ndown`. They become the engine's fermionic quantum
  numbers, conserved or not. Two odd factors must carry the strict
  convention's Jordan–Wigner string `[i, j)` that `ALPS::models` emits; odd
  terms cannot be represented. Under the lattice convention every operator
  is bosonic.
- At most four quantum numbers are tracked: the sector's charges plus the
  fermionic ones. Charges and targets are half-integers in `short` range;
  an empty sector or more targets than the sector holds are rejected before
  the engine starts.
- Real arithmetic needs real coefficients and factor matrices (`Sp Sm` rather
  than `Sy Sy`); complex ones select `Arithmetic::complex`.
- Measured values are real parts, so a non-Hermitian operator such as
  `c†_i c_j` reports the real part of its expectation value. The components
  of one column are measured one by one and summed. A two-site product
  inside one block is measured for a second target from a stored block
  operator, which the solver keeps only when several targets and two-site
  measurements meet.
- Solves are serialized process-wide (the engine keeps global state) and may
  not be nested; each solve owns a private scratch directory under
  `scratch_parent` (default: the system temporary directory) that is removed
  when it returns, and restores the engine's state afterwards, including
  after an exception.

## The `dmrg` application

The XML application `dmrg` (`src/apps/dmrg/dmrg`) reads its parameters,
builds the catalog Hamiltonian and the measured operators with
`alps::models::CatalogModel` (the importer `ALPS::ed` is checked against,
`<alps/model/catalog.hpp>`) and solves them with this component; it is built
only when `ALPS_ENABLE_CLASSIC_DMRG` is on. Parameters and output keep their
names and layout: `SWEEPS`, `STATES` (2 × `SWEEPS` values) or `MAXSTATES`
(ramped) or `NUMSTATES` (fixed), `NUM_WARMUP_STATES`, `TRUNCATION_ERROR`,
`LANCZOS_TOLERANCE`, `NUMBER_EIGENVALUES`, `VERBOSE`, `COMPLEX`,
`TEMP_DIRECTORY`, `DISORDERSEED` and the `MEASURE_*` definitions. A charge is
conserved when `CONSERVED_QUANTUMNUMBERS` names it and `<name>_total` is set;
names the model lacks are ignored. Compared with the application's former
copy of the engine:

- `MEASURE_LOCAL` of a bond operator reports the whole operator, where the
  old code kept only its first term (`Sz Sz`, a third of an isotropic
  exchange).
- `MEASURE_CORRELATIONS` averages over all pairs at each distance on
  periodic lattices, where the old code reported the pairs of site 0
  divided by the multiplicity; on ladders and other lattices whose pair
  order differs from the distance order the values are filed under their
  own labels.
- Two-site values of excited targets are corrected. Energies are equal or
  lower, and the iteration count can differ by one.
- `START_SWEEP` (resuming an interrupted run) is rejected. Scratch files go
  to a private directory below `TEMP_DIRECTORY` that is removed afterwards.
- Lattices of fewer than four sites, more than four tracked charges and
  measurement terms of odd fermion parity are rejected with a message.
- Energies come out equal or lower than the old code's, which stopped its
  Lanczos solves on an energy stall: the solves of the final sweep now
  certify their states ([below](#the-lanczos-stop-rule)). Runs take about as
  long: on an Apple M-series laptop, a spin-1 chain of 64 sites at
  `MAXSTATES=300` takes 17.5 s against 15.8 s for the old code (the certified
  final sweep is the difference), a Hubbard chain of 32 sites at
  `MAXSTATES=400` 6.4 s against 7.5 s, and the dmrg-06 tutorial at 32 sites
  with correlations 3.2 s against 3.1 s.

## The Lanczos stop rule

`Options::lanczos_tolerance` (engine default 1e-7) bounds two quantities of a
solve. Every solve ends once each targeted Ritz energy has changed by at most
the tolerance since the previous Lanczos iteration, as the legacy engine
stopped, and once each targeted residual `|b_k z_k(n)|` is below the square
root of the tolerance, the accuracy a settled energy implies and the guard
against a stalled energy. The solves of the final sweep, whose states are
measured and whose last state's energy is reported, drive every residual
down to the tolerance itself. During the regular sweeps the transformed seed
of the previous iteration is usually accepted after the minimum of three
Lanczos iterations; the certified solves of the final sweep take ten to
seventy at 300 kept states. Certifying every solve (the alpsx rework's rule)
doubled the Lanczos iterations of a run for energies within the tolerance of
these; the legacy stall rule alone left the spin-1 benchmark above 9e-6
above them. A Krylov space exhausted at 99 vectors or broken down ends the
solve as before.

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
disorder state, HDF5 output and run recipe are not ported. Written here and
not taken from alpsx, which fed dmtk from catalog operators: the search for
a set of fermionic carrier charges, factors in either order, the
normalization and combination of terms, the dense Hermiticity check and
the component-wise measurement. The `dmrg` application's own copy of dmtk
was removed when the application moved onto this component.

Engine changes made here, each marked in the affected file: the dense block
products call the platform BLAS through `dmtk/blas.hpp` (Fortran LP64 ABI,
as `ALPS::ed` calls LAPACK) in place of the rework's Eigen views, as the
legacy engine did through the Boost bindings; the Lanczos stop rule
([above](#the-lanczos-stop-rule)); the Lanczos recurrence runs in place and
files only the values of its Krylov vectors (the rework did so for the
reorthogonalized solves); `Vector` and `Matrix` stream their elements in one
call instead of one per scalar (same file layout); and `VectorState::condense`
keeps the merge permutations of the layouts it last saw instead of
rebuilding them for every product.

## Tests

`ctest -L classic_dmrg`: energies against `ALPS::ed` with enough kept states
for exact sweeps (open Heisenberg chain and ladder, spin-1 chain with
anisotropy, sector targets, an unrestricted solve, three targeted states,
hard-core bosons, the Hubbard chain at and away from half filling, frustrated
spinless fermions on a next-nearest-neighbour chain, complex hopping phases,
truncated ramped sweeps converging from above, iteration observation and
abort); measurements against `ALPS::ed` (site and pair spin terms, fermion
densities and hopping correlations, split charge components, scalar terms,
one row per target); equivalent operator forms (reversed factors, scales in
the matrices, projector-form terms, `Nup` and `Ndown` as fermionic charges);
the private scratch directory; and the rejected inputs.
The label also covers XML catalog models solved against `ALPS::ed` (wrap
bonds, spin 1, ladders, spinless fermion rings, Hubbard chains, rings and
ladders, the t-J chain, correlations in both site orders). The SDK consumer
contract solves the open 6-site chain through the installed component.
