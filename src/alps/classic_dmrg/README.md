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

## Molecular Hamiltonians (`<alps/dmrg/molecular.hpp>`)

`solve_molecular(MolecularInput, options, callback)` runs the engine on a
restricted, real ab initio Hamiltonian given by its spatial-orbital integrals
in chemist notation,
`H = Σ h[p,q] c†[p,σ] c[q,σ] + ½ Σ (pq|rs) c†[p,σ] c†[r,τ] c[s,τ] c[q,σ] + core_energy`,
with one spinful-fermion site per orbital in input order and `alpha_electrons`
and `beta_electrons` fixing the sector. It returns the energy (including
`core_energy`), the spin-resolved one-particle reduced density matrices, the
two-particle ones in the PySCF convention
(`dm1[p,q] = <c†[q] c[p]>`, `dm2[p,q,r,s] = Σ <c†[p,σ] c†[r,τ] c[s,τ] c[q,σ]>`),
spin-summed (`two_rdm`) and by spin block (`two_rdm_aa`, `two_rdm_ab` with
σ = α and τ = β, `two_rdm_bb`; `two_rdm` is `aa + ab + ab.transpose(2,3,0,1)
+ bb` to rounding, and QDK's `aaaa`, `aabb` and `bbbb` blocks are these), the
iterations and the final truncation error. The engine works from the
integrals themselves, and the RDMs are read from the final state by a walk
over the chain, not measured as engine terms (both below). The
sweep, truncation, Lanczos, scratch and callback options are those of
`solve()`; `sector` must be empty, `num_eigenvalues` one and arithmetic real.
The integrals must be finite and have the real-orbital symmetries (`h`
symmetric, `(pq|rs)` eightfold) after `integral_cutoff` drops entries of
smaller magnitude. Fewer than four orbitals
are padded with one-dimensional neutral sites, which add no states.

The result keeps no MPS and no checkpoint, so everything is measured from
the final sweep's state and the block rotations and operators it leaves in
the scratch directory. Besides the RDMs, `MolecularInput::observables` names
`models::TermSum`s on `spinful_fermion()` sites (site `i` is orbital `i`)
with real coefficients and at most four sites per product; each comes back
in `MolecularResult::observables` as its expectation value, and products
that change N or Sz are zero without being measured. With
`MolecularInput::ci_vector`, the solver also reconstructs the CI vector of
the final state from the stored block rotations and the final wavefunction:
`ci_alpha`, `ci_beta` (occupation words, bit `p` for orbital `p`) and
`ci_coefficients` list every determinant of the sector, zeros included,
ordered by alpha then beta word, in ALPS mode order (orbitals ascending,
the down mode before the up mode of an orbital), normalized, with the
largest amplitude positive. Its 1-/2-RDMs are those of the run, truncated or
not. The sector grows as `C(n, n_α) C(n, n_β)` (63504 determinants at 10
orbitals and half filling, where the reconstruction adds about a second to
a 30 s run); it is limited to 64 orbitals.

The solver is meant for small and medium active spaces and for checking
other codes. Its engine, `dmtk::MolecularSystem`
(`src/dmtk/molecular_system.hpp`), keeps for every block the normal and
complementary operators of quantum-chemistry DMRG (S. R. White and
R. L. Martin, J. Chem. Phys. 110, 4127 (1999); the operator classes follow
G. K.-L. Chan and M. Head-Gordon, J. Chem. Phys. 116, 4462 (2002)): the
block Hamiltonian, the creation and annihilation operators of the block's
spin orbitals and their pair products, and for every spin orbital outside
the block the pair operators `P` and `Q` and the three-operator sums `S`
that contract the integrals with the block's own operators. A step builds
these for the two enlarged blocks from the stored ones and the new site's,
applies the Hamiltonian as about `2n²` products of one left and one right
operator, and rotates the grown block's operators into the kept states;
the generic engine path rebuilt sums of the `n⁴` integral terms at every
step. With 32 kept states and two sweeps on dense random integrals at half
filling, an Apple M-series laptop shared with other jobs took 0.3 s for 6
orbitals, 1.3 s for 8, 3.7 s for 10, 8.1 s for 12, 17 s for 14 and 32 s for
16, where the generic path took 1.8, 10, 36 and 128 s for 6 to 12 orbitals
and 10.5 minutes for 14 in the same runs; at 64 states 14 orbitals took 30 s
and 16 took 57 s, with a peak of 3.8 GB of memory. Both paths then measured
the RDMs as engine terms.

They no longer do. A measured term stores, in every block the final sweep
grows, the products of its operators that lie inside the block, and for the
RDMs these grow as `n⁴m²`. Even after the final sweep kept only the closing
configuration's left block ([engine changes](#origin-and-provenance)), that
block held 0.76 GB of the 0.9 GB scratch peak of 12 orbitals at 256 states.
Instead, `dmtk::MolecularSystem::measure_products` walks the wavefunction
centre from the final sweep's closing configuration to the left end of the
chain. At each site the left block is the one the final left-to-right half
grew, whose `a`, `a†` and pair operators the engine filed then; the right
block is built on the way, as an orthonormal basis of the state's rows over
the site and the previous right block (Householder QR, nothing truncated, so
the walk carries the final state exactly), together with its annihilators.
Each RDM element is read where the centre sits on its third-lowest orbital,
so the left block holds at most two of its operators and the right block at
most one. The walk keeps one left block's filed operators, the right block's
annihilators and the centre at a time, and costs about one sweep's rotation
of single operators. The caller's observables remain engine terms. Measured
back to back on random integrals at half filling, 12 orbitals at 256 states
and four sweeps went from 99 s to 83 s, from 3.6 GB to 0.88 GB of memory and
from 0.86 GB to 0.18 GB of scratch, and the measurement after the last
iteration from 1.9 s to 0.26 s; 16 orbitals at 64 states and two sweeps
went from 38 s to 19 s, from 2.8 GB to 0.19 GB of memory, from 0.32 GB to
0.04 GB of scratch and from 3.1 s to 0.1 s. What is left is the engine's own block
operators, which grow as `n²m²`: 12 orbitals (a hydrogen chain) at 512
states peak at 2.4 GB of memory and 0.73 GB of scratch with a 0.8 s
measurement, 16 orbitals at 256 states and two sweeps at 1.8 GB and 0.49 GB
(0.8 s), and 10 orbitals with all 1024 states at 3.9 GB and
1.1 GB. alpsx's LiH case (six orbitals, two electrons
per spin) at 128 states and six sweeps took 0.4 s against 2.8 s on the
generic path; alpsx's version of the solver took 3.0 s when the generic
path took 2.3 s. With the walk it takes 0.16 s, against 0.26 s with the
RDMs as engine terms in the same runs.

Converged runs agree with the generic path to rounding. With every state
kept, the energies agree with PySCF's full CI to 2e-14 and the RDMs to 5e-13
(6 and 8 orbitals at half filling, 10 orbitals with two electrons of one
spin and one of the other), and with the generic path's to 6e-15 and 8e-14
(8 and 10 orbitals); while the engine was written, every Hamiltonian product
agreed with the generic path's to 4e-15, from 1 to 10 orbitals and through
warmup, sweeps and final sweep. Runs far from
convergence land elsewhere, as the generic path's own do when the integrals
change at 1e-13. The engine uses the integrals averaged over the eightfold
symmetry that `solve_molecular` checks to 1e-12, which makes the operator
exactly Hermitian.

The RDM walk reads the same state as the measured terms it replaced, so
the two agree to rounding however far the run is from convergence: run side
by side in one process, every element of the 1-RDMs and of the four 2-RDM
spin blocks agreed to 3e-15, from 1 to 10 orbitals, with more alpha than
beta electrons, more beta than alpha or equal numbers, and from 4 to 400
kept states. Against the measured terms of the previous build the RDMs
agree to 2e-14, up to 16 orbitals and 256 states and including the
converged twelve-atom hydrogen chain below, with the spin blocks compared
up to 16 orbitals; the energies and observables are bit for bit unchanged,
and so are the CI vectors (up to 12 orbitals, 853776 determinants). With every state kept and a Lanczos
tolerance of 1e-13, the energies agree with dense full CI to 1e-14 and the
RDMs to 9e-14 (6 and 8 orbitals at half filling, 7 with four and three
electrons). For a hydrogen chain of 10 atoms (STO-3G, 1 Å, all 1024 states)
the RDMs agree with PySCF's full CI to 3e-10, the accuracy of that reference
vector (its residual is 6e-10). Twelve atoms in Löwdin-orthogonalized
atomic orbitals at 256 states and six sweeps come within 1.3e-10 of the
full-CI energy, with the RDMs within 7e-7. In the canonical RHF orbitals,
which are delocalized, the same chain at 512 states stays 3.6e-4 above it.

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
  `MAXSTATES=300` takes 16.6 s against 15.8 s for the old code (the certified
  final sweep is the difference), a Hubbard chain of 32 sites at
  `MAXSTATES=400` 6.1 s against 7.5 s, and the dmrg-06 tutorial at 32 sites
  with correlations 2.9 s against 3.1 s.

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

The molecular solver (`include/alps/dmrg/molecular.hpp`,
`src/molecular.ipp`, included at the end of `src/classic.cpp`) is ported from
alpsx branch `WIP/classic-dmrg-chemistry`, commit d4384a7f8; its files name
the adaptations, which concern the headers, the callback and the charge
registration. Its tests are adapted from the same branch
(`tests/dmrg/classic_molecular.cpp`) and from alpsx commit 95cf817d6
(`tests/dmrg/classic_multisite_fermion_lock.cpp`). Written here, not in
alpsx: the spin blocks of the two-particle RDM, the observables and the CI
reconstruction, and their tests.

Engine changes made here, each marked in the affected file: the dense block
products call the platform BLAS through `dmtk/blas.hpp` (Fortran LP64 ABI,
as `ALPS::ed` calls LAPACK) in place of the rework's Eigen views, as the
legacy engine did through the Boost bindings; the Lanczos stop rule
([above](#the-lanczos-stop-rule)); the Lanczos recurrence runs in place and
files only the values of its Krylov vectors (the rework did so for the
reorthogonalized solves); `Vector` and `Matrix` stream their elements in one
call instead of one per scalar (same file layout); `VectorState::condense`
keeps the merge permutations of the layouts it last saw instead of
rebuilding them for every product; and the per-step bookkeeping of
`System` (`build_composite_operators`, `hint`, `rotate_terms`) finds block
operators through an index built once per call instead of a list scan per
lookup, sorts a term's factors into blocks before copying any of them, and
groups the terms whose block operators are summed into one composite
operator through a table instead of comparing every pair of terms, which
was quadratic in the number of terms. The grouping keeps the summed
operator but can sum it in a different order (the replaced scan stepped
over the term after each one it merged): converged runs agree with the old
code to rounding, and a run far from convergence can land elsewhere, as
the old code's own does when the integrals change at 1e-13 (at 8 and 10
orbitals with 32 states and two sweeps, by 1e-3 and 1e-2 in the energy).
The molecular solver adds its measured terms (now only the caller's
observables) for the final sweep only, instead of rotating them through
every step. With stored products and a full sweep
(the molecular solver, and the `dmrg` application with several targets and
two-site measurements), `final_sweep` builds measurement operators only
where the closing measurement reads them. It skips the right blocks of its
right-to-left half and the left block grown at its last step, and it
removes each left block's file once the next step has read it; `read_block`
then refuses that block. Measured values are bit for bit the same: the
molecular energies and RDMs at 10, 12 and 16 orbitals (the RDMs were then
measured terms), and the `dmrg`
application's XML output for a 12-site spin chain with two targets and
correlations.
`System` takes its Hamiltonian through three virtual hooks,
`init_hamiltonian` (the end of `init_iteration`), `hamiltonian_product` (the
product of the Lanczos and one-state solves) and `rotate_hamiltonian` (in
`rotate`), whose defaults are the code that stood there; the molecular
solver's `dmtk::MolecularSystem` overrides them and is written here, as is
its RDM walk `measure_products`, not taken from alpsx or the legacy engine. `rotate_corr` and `measure` also find
block operators through the per-call index, and `rotate_corr` answers
whether the new block already holds an operator from a key set instead of a
scan of its list, which was quadratic in the number of measured terms.

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
The molecular solver (`classic_dmrg_molecular`, also labelled `chemistry`)
is checked against an occupation-number exact diagonalization for one to five
and for seven orbitals, with and without padding, energies and every RDM element; against
the same integrals written as an `ALPS::models` Hamiltonian and solved by
`ALPS::ed`; the spin blocks, observables of one to four sites and the CI
vector, sign by sign, against the same diagonalization for one to six
orbitals; under truncation (variational energy, RDM traces, the CI vector's
norm, and its 1-RDMs and 2-RDM spin blocks against the run's); for callback
observation and abort with the scratch directory removed; and for the rejected
inputs. `classic_dmrg_multisite` locks the engine's long-range hopping and
four-site fermion products, which the molecular solver used before it had
its own operator engine, against an exact diagonalization for 4, 6 and 8
sites.
The label also covers XML catalog models solved against `ALPS::ed` (wrap
bonds, spin 1, ladders, spinless fermion rings, Hubbard chains, rings and
ladders, the t-J chain, correlations in both site orders). The SDK consumer
contract solves the open 6-site chain through the installed component.
