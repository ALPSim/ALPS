# `ALPS::classic_dmrg` — classic DMRG on compiled models

Finite-system, two-site density-matrix renormalization group for models built
with [`ALPS::models`](../models/README.md): Adrian Feiguin's dmtk engine
driven by the numeric local terms of a `CompiledOperator` on its
`PhysicalBasis`. The XML `dmrg` application runs on this component
([below](#the-dmrg-application)). The engine warms up, sweeps with a
kept-state schedule, targets one or several lowest states of one conserved
sector and measures compiled operators on them, or evolves the ground state
in real or imaginary time ([below](#time-evolution)). The component depends on
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
| `Options` | `sweeps` (`fixed_sweeps`, `ramped_sweeps` or an explicit `SweepStates` list; the right-to-left half sweep runs first), `warmup_states`, `warmup_configuration` ([below](#unevenly-distributed-ground-states)), `num_eigenvalues`, optional `truncation_error` with `truncation_max_states`, optional `lanczos_tolerance` ([below](#the-lanczos-stop-rule)), `arithmetic` (real or complex), `sector` (name → total of a declared conserved charge), `scratch_parent`, `log`, `verbose`, `keep_mps` ([below](#matrix-product-states)) |
| `Measurements` | `expectation`: name → operator, one value `<ψ\|O\|ψ>`; `terms`: name → operator, one value per term (one column per site, pair, …) |
| `site_terms(model, op)`, `pair_terms(model, a, b)` | builders over a model's local operator names: one term per site, one term per ordered site pair (i-major, on-site product at i = j); no two columns are ever combined |
| `solve(op, basis, options, measurements, callback)` | the engine on a compiled operator and its physical basis |
| `solve(model, …)` | validates Hermiticity, compiles the `Hamiltonian` under its declared charges and convention, solves on `physical_basis(model)` |
| `Result` | ascending `energies`, `measurements` (rows follow the targets, columns the terms with their factor `sites`), the `iterations` observed (direction, index, energy, truncation error, entropy), `final_truncation_error`, the resolved `options`, with `keep_mps` the `mps` of every target |
| `MatrixProductState`, `MpsTensor` | a targeted state's site tensors (`[left bond][local state][right bond]`, row-major, `real` or `complex` entries) and its `center` |
| `evolve(initial, final, basis, evolution, options, measurements, callback)`, `evolve(initial_model, final_model, …)` | time evolution of the ground state of `initial` under `final` ([below](#time-evolution)) |
| `Evolution` | `time_step`, `steps`, `imaginary`, `states`, optional `truncation_error`, `measure_steps`, `drives`, `keep_mps`, `method` (`EvolutionMethod::trotter` or `targeting`, [below](#time-step-targeting)), `half_sweeps` |
| `Drive` | a time-dependent term: `coefficients[n] * op` during step `n` |
| `EvolutionResult` | `initial` (the preparing solve's `Result`), the measured `steps`, their `times`, `energies` (`<H_final>`) and `measurements` (one row per measured step), the `truncation_errors` of every step, with `keep_mps` the `mps` of every measured step |
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

## Unevenly distributed ground states

The warm-up grows the chain from the left, one site at a time, with a
one-site environment on the right. Each partial chain of `k` sites targets
the share `k / L` of the sector, rounded to the nearest reachable charge.
When the ground state distributes its charges unevenly (a domain wall,
random fields, a trap), those uniform shares misplace them.
The sweeps then move charge only across the cut between the two blocks,
so they correct a misplaced charge only where an exchange across the cut
lowers the energy. With fields alone nothing moves, and the solve keeps an
excited product state. Weak hopping moves charges slowly: with 1e-3 of the
exchange, 24 sweeps still leave three spin pairs misplaced on the wall
below (-18.006 instead of -24). Nothing warns about either case. Measured on 48 sites with fields `-1`
left of site 12 and `+1` from it on (`Sz = -12`), four sweeps of 32 states:

| initial Hamiltonian | default warm-up | with `warmup_configuration` | exact |
| --- | --- | --- | --- |
| fields alone | -8 | -24 | -24 |
| fields and XX exchange `J = 0.1` | -18.557 (-24.00125 after 12 sweeps) | -24.0012508 | -24.0012508 (free fermions) |

`Options::warmup_configuration` names a product configuration instead,
as one local-state index per site, indexing that site's states in the
physical basis. For spin 1/2, state 0 is up and state 1 is down. The
partial chain of `k` sites then targets the charges of the
configuration's first `k` sites. This uses the engine's own custom-target
warm-up (`System::warmup_loop` with a target per length), so no engine
file changes. The configuration only places the charges: the sweeps still
optimize the state, and any configuration with the right charge profile
works. It needs one state per site, and its conserved charges must add up
to the sector; `solve()` and `evolve()` reject anything else. The
`evolve()` preparation uses it too. `ResolvedOptions` reports it.

Without a configuration the default stays the uniform warm-up. The
alternatives the engine already offers each fail elsewhere:

- growing from both ends (`set_grow_symmetric`) finds the centred wall
  but not the one at site 12;
- the warm-up over all charge subspaces (`set_target_subspaces`, random
  states in every subspace) finds both walls but loses random fields;
- random Lanczos starts in every sweep still stop at `E = -14` for the
  48-site centred wall.

A weak coupling added to fields alone does not fix it either: 1e-3 of
the exchange gives -14.009 instead of -24 on the centred 48-site wall
with four sweeps.

## Matrix product states

With `Options::keep_mps`, `Result::mps` holds every targeted state as a
`MatrixProductState`, in energy order. Nothing is approximated beyond the
DMRG itself: finite-system DMRG builds each block from the previous one and
a site through a stored rotation, so the rotations of the final sweep are the
tensors up to its closing configuration, whose two-site wave function an SVD
splits into the remaining two (singular values below 1e-14 of the largest are
dropped). They are read from the scratch files and the engine's wave
functions before the solve returns. The state is mixed canonical with its
norm, one, at `center`: the middle (`L/2`) when the final sweep ends there,
and `L-2` when the engine stores products and sweeps the whole chain
(several targets with two-site measurements). The targets share every tensor
but the two split ones. Local states follow each site's `PhysicalBasis`
order, mapped from the engine's charge-sorted bases; with fermions a product
state is the product of creators in site order applied to the vacuum, the
convention of `ALPS::ed`, so the MPS contracted to a vector is the `ALPS::ed`
eigenvector up to a phase. `solve_molecular` honours `keep_mps` too and
returns the state on the orbitals (padding sites folded into the last one).
pyalps returns the tensors as NumPy arrays, with `to_dense`, `amplitude`,
`overlap` and `norm`. Measurements after the solve are not part of the
engine; a time evolution keeps its states through `Evolution::keep_mps`
([below](#time-evolution)).

## Time evolution

`evolve(initial, final, basis, evolution, options, measurements, callback)`
prepares the ground state of `initial` in `options.sector` by the usual
warm-up and sweeps (the last sweep certifies its Lanczos residuals) and
evolves it under `final` by adaptive time-dependent DMRG (S. R. White and
A. E. Feiguin, Phys. Rev. Lett. 93, 076401 (2004)). This section describes
the default method, `EvolutionMethod::trotter`; time-step targeting
([below](#time-step-targeting)) is the other. Every step is the
second-order Trotter product `exp(-dt A/2) exp(-dt B) exp(-dt A/2)`, with
`A` the bonds `(b, b + 1)` of even `b` and `B` the odd ones, and `dt`
imaginary in real time. Each factor is one half sweep over the whole chain
that applies the factor's bond gates where both sites of a bond are
explicit, the central bond and, while a block is a single site, the edge
bond, so every bond is visited once. The state is carried from step to
step by White's wave-function transformation, no Lanczos solve runs, and
the density matrix of the evolved state chooses the kept states, so the
basis follows the state. Two adjacent `A` half steps are merged into one
half sweep of their gate products, which is exact (drives included); a
measured step closes its `A` half step and is measured in a half sweep of
identity gates (and step 0 in two, which first bring the measured
operators into the blocks the preparation left). A step therefore costs
two half sweeps, a measured one four.

- The gates are built from the compiled terms: a bond's matrix is its
  two-site terms plus the one-site terms of its two sites, each site's
  divided evenly between its bonds (an end site's go wholly to its one
  bond). A product of two odd factors carries the parity of the bond's
  left site in the input state, the strict convention's string as
  `ALPS::ed` applies it. The exponential is taken with Eigen on every
  connected block of the bond matrix, so the gates conserve every charge
  exactly. `final` and every drive must be Hermitian, conserve the
  sector's charges and, for the Trotter method, couple only sites `i` and
  `i + 1`; a coupling of other sites (a periodic chain's closing bond,
  next-nearest neighbours) is rejected with a message naming the sites.
  Time-step targeting takes them.
- Real time runs in complex arithmetic; imaginary time (`imaginary`) runs in
  `options.arithmetic`, applies `exp(-tau H)` and renormalizes after every
  gate, approaching the lowest state of `final` that the initial state
  overlaps, `O(dt^2)` above it in energy (the fixed point of the
  Trotterized propagator).
- `states` caps every truncation of the evolution; with `truncation_error`
  a truncation keeps as many states as that discarded weight needs, up to
  the cap (dmtk adds a tenth). `options` controls the preparation only;
  `num_eigenvalues` must be one, and `Options::keep_mps` is replaced by
  `Evolution::keep_mps`.
- `measure_steps` lists the measured steps, 0 being the initial state;
  empty measures every step. Every measured step reports its time
  (`steps * time_step`), `<H_final>` without the drives, the requested
  measurements as `solve()` takes them, and with `keep_mps` the state as a
  `MatrixProductState` whose centre is where the measuring half sweep
  ended (site 2 or `L - 2`).
- `truncation_errors[n]` is the weight discarded in the half sweeps of step
  `n` (the shared half sweep of two merged half steps counts for the
  earlier step; entry 0 covers the initial measurement); the squared norm
  the state has lost is about their sum.
- A `Drive` adds `coefficients[n] * op` to `final` during step `n`; pyalps
  samples a function of time at the step midpoints `(n + 1/2) dt`, which
  keeps the scheme second order.
- The callback observes the preparation's iterations and then every
  iteration of the evolution, which carry no energy (NaN) but their
  truncation error and entropy.

A product initial state, such as a domain wall prepared as the ground
state of fields alone, needs `Options::warmup_configuration`. Without it,
the preparation can return an excited product state
([above](#unevenly-distributed-ground-states)). The domain-wall tutorial
passes the wall as its configuration.

The tests below compare every step with the dense Trotter product and the
exact evolution on the ED sector basis; with enough kept states the states
agree with the Trotter product to 1e-9. The XX domain wall of 24 sites
follows the exact Bessel solution to 5e-5 up to `t = 4` at `dt = 0.05`,
120 steps with up to 40 states in about 3 s on an Apple M-series laptop.
pyalps binds it as `classic_dmrg.evolve`, with `overlaps` and
`loschmidt_echo`; [tutorials/time-evolution](../../../tutorials/time-evolution/README.md)
has four worked examples.

### Time-step targeting

`Evolution::method = EvolutionMethod::targeting` (pyalps
`method="targeting"`) evolves by time-step targeting (A. E. Feiguin and
S. R. White, Phys. Rev. B 72, 020404(R) (2005)). It has neither the
nearest-neighbour restriction nor a time-step error. Instead of gates, it
applies `final` itself to the superblock state:

- `final` replaces the system's Hamiltonian under its own name, and the
  hooks that build the ground state's block operators build its block
  operators. Two half sweeps do this, one side each, before the first step;
  the first of them targets the state alone. `final` may couple any two
  sites, as the initial Hamiltonian may (next-nearest neighbours, a
  periodic chain).
- At every configuration four Runge-Kutta vectors of `c H`, with
  `c = -i dt` (or `-dt` in imaginary time), give the states at `t + dt/3`,
  `t + 2dt/3` and `t + dt` (the paper's formulas, in `src/tdmrg.hpp`). The
  density matrix targets them and the state with weights 1/3, 1/6, 1/6 and
  1/3, so the kept basis holds the whole interval of the step.
- Each step is `half_sweeps` half sweeps (default 2). At the last
  configuration of the last one, a Krylov exponential of `final` (Lanczos
  with full reorthogonalization, up to 40 vectors) advances the state.
  With enough states this is exact, so a large `dt` costs only the Krylov
  vectors. Imaginary time renormalizes after the step.
- A measured step is measured in one more half sweep. That half sweep
  keeps the four targets and does not advance. Step 0 is measured in the
  second of the two building half sweeps.
- Drives are rejected: the block Hamiltonians hold `final`'s coefficients,
  built once. `truncation_errors[n]` is the squared norm the state itself
  lost in step `n` (the four targets share the kept states, so a
  truncation keeps more states than the state alone would need). The
  iterations report `<psi|H_final|psi>` at their configuration, except the
  first building half sweep's, which report NaN.

Each configuration costs four products with `H`, and the last one of a step
also costs the Krylov vectors. The mixed density matrix keeps more states
than one state needs. The 24-site domain wall of the tutorial, with up to 40 states
and a truncation error of 1e-10, follows the Bessel solution up to `t = 4`
as follows (Apple M-series laptop):

| method | `dt` | steps | time | max deviation |
| --- | --- | --- | --- | --- |
| Trotter | 0.05 | 120 | 2.4 s | 4.8e-5 |
| targeting | 0.05 | 120 | 3.6 s | 2.1e-5 |
| targeting | 0.2 | 30 | 1.4 s | 3.2e-5 |

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
not. The sector grows as `C(n, n_α) C(n, n_β)`: at half filling with 32
kept states and two sweeps, the reconstruction adds 0.03 s to a 1.4 s run
at 10 orbitals (63504 determinants) and 0.35 s to a 2.5 s run at 12
(853776). It is limited to 64 orbitals.

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
`TEMP_DIRECTORY`, `DISORDERSEED`, `SAVE_MPS` and the `MEASURE_*` definitions. A charge is
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
- `SAVE_MPS=1` writes the ground state's `MatrixProductState` to `/mps` of
  the HDF5 result, in the layout of master's
  [ALPSim/ALPS#216](https://redirect.github.com/ALPSim/ALPS/pull/216), which
  `pyalps.loadMPS` reads: `length`, `center`, `tensors/<i>/{shape,values}`,
  `site_type` and the local quantum numbers under `site_basis`. `center` is
  where the final sweep ended (above), not always `L/2`.
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

## Citing

The citation catalog ([CITATIONS.yaml](../../../CITATIONS.yaml), rendered in
[CITATION.md](../../../CITATION.md)) has two components for this library:
`classic_dmrg` for ground states (`solve()`, with the same references as the
`dmrg` application) and `classic_dmrg_evolution` for `evolve()`.
`alps::print_citation_details(out, "classic_dmrg")` (`<alps/utility/citations.hpp>`) prints the complete
references.

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
rebuilding them for every product; `System::init_iteration` (through its
default `init_hamiltonian`) builds the interaction and composite terms once
per step instead of twice (a guarded
copy of the same calls ran first on the same inputs, as in the legacy
engine); and the per-step bookkeeping of `System` (`build_composite_operators`,
`hint`, `rotate_terms`) finds block
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
its RDM walk `measure_products`, not taken from alpsx or the legacy engine.
`System::signal_emit(std::string)` returns `true` for an unknown signal
instead of falling off its end, the fix of
[ALPSim/ALPS#232](https://redirect.github.com/ALPSim/ALPS/pull/232) to the
legacy engine.

The time evolution is written here: the bond matrices, the schedule and
the result in `src/classic.cpp`, and the engine side, `detail::BondGate`,
the gate exponentials and their application to the engine's
charge-sorted wave functions, and `detail::EvolvingSystem`, in
`src/tdmrg.hpp`. It is the method of White and Feiguin (PRL 93, 076401);
the sweep schedule (even bonds `dt/2`, odd bonds `dt`, even bonds `dt/2`,
a measuring half sweep of identity gates; gates on the central bond and on
an edge bond whose block is one site; no Lanczos solve) was read from
Adrian Feiguin's teaching program
[afeiguin/dmrg-tutorials](https://github.com/afeiguin/dmrg-tutorials)
(`06-heischain_time/heischain_time.cc`), a dense spin-1/2 Heisenberg code
without quantum numbers. That repository has no license, so it served as a
description of the algorithm only and no code was copied or translated
from it; maintainers may want to confirm that this use needs no
permission. Time-step targeting in the same file was written from the
paper (Feiguin and White, PRB 72, 020404(R)), with its Runge-Kutta targets
and weights; the Krylov step is a textbook Lanczos exponential. No dmtk
file changed: `EvolvingSystem` derives from `dmtk::System` and runs its
warm-up and sweeps. For the Trotter method it disables the two Hamiltonian
hooks while the state evolves, so only the measured operators are built
and rotated. For time-step targeting it keeps them and replaces the
system's Hamiltonian by the final one.

The MPS export (`src/mps.hpp`) is adapted from master's SAVE_MPS
([ALPSim/ALPS#216](https://redirect.github.com/ALPSim/ALPS/pull/216), Marcus
Rosales, commits 46c76b1d2 and 35cea33e7): `DMRGTask::extract_mps` and
`split_two_site` of the legacy application, which assembled its ground state
from the same rotations. Here it reads the closing configuration and every
target from the engine in memory instead of the ground-state file at `L/2`,
splits with Eigen's SVD instead of LAPACK, and maps local states through the
engine's sorted bases as the molecular CI walk does. The application's
`/mps` output keeps #216's layout, and `pyalps.loadMPS` is #216's. `rotate_corr` and `measure` also find
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
the private scratch directory; the rejected inputs; and the kept MPS
(`classic_dmrg_mps`), contracted to a vector and compared with the
`ALPS::ed` eigenvector of each target for spin-1/2 and spin-1 chains, an odd
length, four sites, an unrestricted solve, the Hubbard chain, frustrated
spinless fermions, complex hopping phases, excited targets and stored
products, with its canonical form checked, also under truncation.
The time evolution (`classic_dmrg_evolution`) compares every measured step
with the dense Trotter product and the exact `exp(-iHt)` on the `ALPS::ed`
sector basis: an XXZ quench in a field, odd and four-site chains, the XX
domain wall, a Hubbard quench, interacting spinless fermions and complex
hopping phases (the gates' parity strings and phases), with energies,
site and pair measurements, the kept states and the Loschmidt amplitude;
it checks an identity and an eigenstate evolution, the `dt^2` shrinking of
the Trotter error, energy conservation, imaginary time in real and complex
arithmetic (against the Trotter product and the ED ground energy), drives,
the measured-step bookkeeping, the reported truncation and the rejections.
Its time-step targeting tests compare every step with `exp(-iHt)` to 1e-9,
and every iteration's energy with the conserved `<H_final>`:
- fields to the XX chain, whose block operators the evolution builds;
- next-nearest-neighbour and periodic final Hamiltonians;
- a Hubbard quench;
- complex hopping;
- four sites;
- imaginary time in real and complex arithmetic.

Two further tests check the method's distinguishing properties. At
`dt = 0.4`, targeting is exact where the Trotter product misses by more
than 1e-2. A truncated product-state domain wall stays within 1e-4 of the
exact state, where targeting the state alone misses by 0.025 after one
step.
`classic_dmrg_warmup` solves 48-site chains whose uniform warm-up
misplaces charges with the default four sweeps, each failing without its
warm-up configuration: field-only walls at the centre and at a quarter
(energy and every `Sz`), irregular fields (against the sorted optimum),
the XX chain in a wall and in irregular fields (against free fermions),
a Hubbard potential step with two charges, and an `evolve()` preparation.
The molecular solver (`classic_dmrg_molecular`, also labelled `chemistry`)
is checked against an occupation-number exact diagonalization for one to five
and for seven orbitals, with and without padding, energies and every RDM element; against
the same integrals written as an `ALPS::models` Hamiltonian and solved by
`ALPS::ed`; the spin blocks, observables of one to four sites and the CI
vector, sign by sign, against the same diagonalization for one to six
orbitals; under truncation (variational energy, RDM traces, the CI vector's
norm, and its 1-RDMs and 2-RDM spin blocks against the run's); for callback
observation and abort with the scratch directory removed; and for the rejected
inputs; and its kept MPS against the CI vector of the same run on every
determinant, padded, exact and truncated. `classic_dmrg_multisite` locks the engine's long-range hopping and
four-site fermion products, which the molecular solver used before it had
its own operator engine, against an exact diagonalization for 4, 6 and 8
sites. With the `dmrg` application built, `ctest -L dmrg` also runs
`dmrg_save_mps`, which reads back the `/mps` of a spin chain and of a
Hubbard chain with two targets and correlations: shapes, canonical form,
norm, the charge total of every product state with weight, and the centre.
The label also covers XML catalog models solved against `ALPS::ed` (wrap
bonds, spin 1, ladders, spinless fermion rings, Hubbard chains, rings and
ladders, the t-J chain, correlations in both site orders). The SDK consumer
contract solves the open 6-site chain through the installed component.
