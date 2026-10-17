// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/include/alps/dmrg/classic.hpp. Adaptation: the
// solver takes an ALPS::models CompiledOperator with its PhysicalBasis (or a
// Hamiltonian) instead of a catalog model; measurements are compiled
// operators instead of catalog expressions; the catalog snapshot, random
// disorder state and run recipe are not ported.
//
// Classic (two-site, finite-system) DMRG on compiled models: Adrian Feiguin's
// dmtk engine, driven by the numeric local terms of ALPS::models. The engine
// runs a warm-up followed by finite sweeps, targets one or several lowest
// states of one conserved sector and measures compiled operators on them.
#pragma once

#include <alps/classic_dmrg_export.h>
#include <alps/models/compiled.hpp>
#include <alps/models/hamiltonian.hpp>
#include <alps/models/physical_basis.hpp>

#include <complex>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <iosfwd>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace alps::dmrg::classic {

enum class Arithmetic { real, complex };
enum class Direction { left_to_right, right_to_left };

// Kept states of one sweep: the right-to-left half sweep runs first.
struct SweepStates {
    std::size_t right_to_left;
    std::size_t left_to_right;
    bool operator==(SweepStates const& o) const {
        return right_to_left == o.right_to_left && left_to_right == o.left_to_right;
    }
    bool operator!=(SweepStates const& o) const { return !(*this == o); }
};

ALPS_CLASSIC_DMRG_DECL std::vector<SweepStates> fixed_sweeps(std::size_t sweeps, std::size_t states);
// The historical MAXSTATES integer ramp: floor((half_sweep + 1) * maximum / (2 * sweeps)).
// A ramp that would retain zero states is rejected.
ALPS_CLASSIC_DMRG_DECL std::vector<SweepStates> ramped_sweeps(std::size_t sweeps, std::size_t maximum);

struct Options {
    std::vector<SweepStates> sweeps = fixed_sweeps(4, 40);
    std::size_t warmup_states = 20;
    // A product configuration, one local-state index per site (into that
    // site's states in the physical basis), whose prefix charges the warm-up
    // targets: its growing chain of k sites takes the charges of the
    // configuration's first k sites instead of the uniform share of the
    // sector. Give it when the charges of the ground state are distributed
    // unevenly (domain walls, random fields, traps). The
    // uniform warm-up misplaces them there, and the sweeps move charge only
    // across the cut: a Hamiltonian of fields alone keeps the misplaced
    // product state, and weak hopping corrects it slowly.
    // Its conserved charges must add up to the sector. Empty: uniform shares.
    std::vector<std::size_t> warmup_configuration;
    // Lowest states targeted together; the ground state is always the first.
    std::size_t num_eigenvalues = 1;
    // Keep as many states as this truncation error needs, up to
    // truncation_max_states (which requires truncation_error), instead of the
    // scheduled counts.
    std::optional<double> truncation_error;
    std::optional<std::size_t> truncation_max_states;
    // Lanczos stop rule (engine default 1e-7): a solve ends once every
    // targeted Ritz energy has changed by at most this between iterations
    // and every targeted residual is below its square root; the solves of
    // the final sweep, whose states are measured, drive every residual down
    // to the tolerance itself.
    std::optional<double> lanczos_tolerance;
    Arithmetic arithmetic = Arithmetic::real;
    // Fixed totals of declared conserved charges (name -> value); every
    // charge of the sector must be declared by the physical basis. Omitted
    // charges remain mixed.
    std::map<std::string, double> sector;
    // Parent of the private scratch directory; empty selects the system
    // temporary directory. The directory is removed when the solve ends.
    std::filesystem::path scratch_parent;
    // Engine diagnostics; nullptr discards them. The caller keeps ownership.
    std::ostream* log = nullptr;
    int verbose = 0;
    // Also return every targeted state as a matrix product state (Result::mps).
    bool keep_mps = false;
};

// One site tensor of a matrix product state: entry (l, s, r) is at
// (l * local + s) * right + r, for left bond l, local state s in the site's
// PhysicalBasis order and right bond r. The outer bonds of the chain have
// dimension one. `real` holds the entries in real arithmetic, `complex` in
// complex arithmetic; the other is empty.
struct MpsTensor {
    std::size_t left = 0;
    std::size_t local = 0;
    std::size_t right = 0;
    std::vector<double> real;
    std::vector<std::complex<double>> complex;
};

// A targeted state as the engine holds it after the final sweep, which is
// no approximation of it: the block rotations of that sweep are the tensors
// up to its closing configuration, and the closing two-site wave function is
// split by an SVD, which drops singular values below 1e-14 of the largest.
// Mixed canonical with the norm, one, at `center`: the tensors before it are
// left-orthonormal, those after it right-orthonormal. The centre is where the
// final sweep ends, at site L/2, or at site L-2 when the engine stores
// products (several targets with two-site measurements); a time evolution's
// states (EvolutionResult::mps) are taken where its measuring sweep ends, at
// site 2 or site L-2. The amplitude of a
// product state is the product of its tensor slices; with fermions, the
// product state is the product of creators in site order applied to the
// vacuum, as in ALPS::ed.
struct MatrixProductState {
    std::vector<MpsTensor> sites;
    std::size_t center = 0;
};

// Operators measured on every targeted state. `expectation` reports one value
// per operator, <psi|O|psi>; `terms` reports the expectation value of every
// term of the operator separately (one column per term, in term order), for
// example one column per site or per site pair. Every term must act on at
// most two sites and have even fermion parity. Values are real parts, so a
// non-Hermitian operator such as c+_i c_j reports Re <psi|O|psi>.
struct Measurements {
    std::map<std::string, models::CompiledOperator> expectation;
    std::map<std::string, models::CompiledOperator> terms;
};

enum class MeasurementKind { expectation, terms };
struct Measurement {
    MeasurementKind kind;
    // One entry per column: the factor sites of the measured term (empty for
    // an expectation value or a scalar term).
    std::vector<std::vector<int>> sites;
    // Rows follow the targeted states in energy order; columns follow `sites`.
    // The engine reports real parts, including in complex arithmetic.
    std::vector<std::vector<double>> values;
};

struct Iteration {
    Direction direction;
    std::size_t index;
    double energy;
    double truncation_error;
    double entropy;
};

// The scientific controls that took effect, without log stream, scratch path
// or callback. The Lanczos tolerance is read back from the engine.
struct ResolvedOptions {
    std::vector<SweepStates> sweeps;
    std::size_t warmup_states = 0;
    std::vector<std::size_t> warmup_configuration;
    std::size_t num_eigenvalues = 0;
    std::optional<double> truncation_error;
    std::optional<std::size_t> truncation_max_states;
    double lanczos_tolerance = 0.;
    Arithmetic arithmetic = Arithmetic::real;
    std::map<std::string, double> sector;
};

struct Result {
    // Ascending energies of the targeted states.
    std::vector<double> energies;
    std::map<std::string, Measurement> measurements;
    std::vector<Iteration> iterations;
    double final_truncation_error = 0.;
    ResolvedOptions options;
    // With Options::keep_mps: one per targeted state, in energy order. The
    // targets share the tensors outside the two sites split at the closing
    // configuration.
    std::vector<MatrixProductState> mps;
};

// Observes every completed DMRG iteration and may throw to abort the solve.
// Solves are serialized process-wide; the callback must not start another one.
using IterationCallback = std::function<void(Iteration const&)>;

// Solve a compiled Hermitian operator on its physical basis. The basis must
// list every site of the operator with its charges and parities; a fermionic
// basis (any odd-parity state under the strict convention) must carry integer
// charges whose summed parity is the fermion parity on every site: the
// particle number N, or else the smallest such set (the XML catalog's Nup and
// Ndown), which become the engine's fermionic quantum numbers. Terms are
// products of at most two local factors on distinct sites, in either site
// order, with the Jordan-Wigner strings ALPS::models emits; like terms are
// combined (also terms that differ only by a scale moved between coefficient
// and factor matrices) and the sum must be Hermitian, which is checked on the
// dense site-pair matrices. At most four quantum numbers (the sector's charges
// plus the fermionic ones) can be tracked. At least four sites are needed.
// Real arithmetic requires real coefficients and factors.
ALPS_CLASSIC_DMRG_DECL Result solve(models::CompiledOperator const& hamiltonian,
    models::PhysicalBasis const& basis, Options const& options = {},
    Measurements const& measurements = {}, IterationCallback callback = {});

// Validate Hermiticity, compile the model under its declared charges and
// fermion convention, and solve it on its physical basis.
ALPS_CLASSIC_DMRG_DECL Result solve(models::Hamiltonian const& model, Options const& options = {},
    Measurements const& measurements = {}, IterationCallback callback = {});

// ------------------------------------------------------------ time evolution
//
// Adaptive time-dependent DMRG (S. R. White and A. E. Feiguin, Phys. Rev.
// Lett. 93, 076401 (2004)) on the same engine: the initial state is the DMRG
// ground state of an initial Hamiltonian in a sector, and it evolves under a
// final Hamiltonian on the same sites by a second-order Trotter
// decomposition into nearest-neighbour bond gates, exp(-i dt h_b) in real
// time or exp(-dt h_b) in imaginary time, each applied where both sites of
// its bond are explicit in the sweeping configuration while the block bases
// adapt to the evolved state (no Lanczos solve). One step is the gates of the even bonds for
// dt/2, of the odd bonds for dt and of the even bonds for dt/2, each set in
// one half sweep; the closing and opening even half steps of consecutive
// steps share a half sweep unless the state between them is measured, which
// takes one more half sweep with identity gates. h_b is the bond's two-site
// terms plus each site's one-site terms divided among its bonds.

// A time-dependent term: coefficients[n] * op is added to the final
// Hamiltonian during step n (from time n dt to (n + 1) dt). Sampling a
// protocol f at the step midpoints, f((n + 1/2) dt), keeps the scheme second
// order. The operator follows the final Hamiltonian's rules.
struct Drive {
    models::CompiledOperator op;
    std::vector<double> coefficients;
};

struct Evolution {
    double time_step = 0.05;
    std::size_t steps = 0;
    // exp(-tau H), the state renormalized after every gate: it approaches
    // the lowest state of the final Hamiltonian in the initial state's
    // sector that the initial state overlaps.
    bool imaginary = false;
    // The most states any truncation keeps; with truncation_error, as many
    // as that discarded weight needs (dmtk adds a tenth), up to this.
    std::size_t states = 100;
    std::optional<double> truncation_error;
    // Steps after which the state is measured, 0 being the initial state;
    // empty measures every step from 0 to `steps`.
    std::vector<std::size_t> measure_steps;
    std::vector<Drive> drives;
    // Also return the state at every measured step (EvolutionResult::mps).
    bool keep_mps = false;
};

struct EvolutionResult {
    // The DMRG of the initial Hamiltonian that prepared the initial state:
    // its energy, iterations and resolved options (with the arithmetic the
    // whole run used). Its measurements and mps are empty: measure step 0
    // for the initial state's.
    Result initial;
    // The measured steps, ascending, and their times (steps * time_step;
    // the imaginary time tau in imaginary time).
    std::vector<std::size_t> steps;
    std::vector<double> times;
    // <psi|H_final|psi> at the measured steps, without the drives.
    std::vector<double> energies;
    // One row per measured step (Measurement::values[row][column]).
    std::map<std::string, Measurement> measurements;
    // The weight discarded in preparing each step's state, summed over its
    // truncations: entry 0 covers the transformation of the initial state
    // and its measurement, entry n the half sweeps that apply step n (the
    // shared half sweep goes to the earlier step) and measure it. The
    // squared norm the state has lost is about their sum.
    std::vector<double> truncation_errors;
    // With Evolution::keep_mps: the state at every measured step.
    std::vector<MatrixProductState> mps;
};

// Evolve the ground state of `initial` (in options.sector) under `final`.
// The final Hamiltonian and every drive must be Hermitian, conserve the
// sector's charges and couple only sites i and i + 1 (chain order); the
// initial Hamiltonian is anything solve() takes. `options` controls the
// initial state's DMRG (its last sweep certifies the Lanczos residuals);
// real-time evolution always runs in complex arithmetic, imaginary time in
// options.arithmetic. num_eigenvalues must be one and keep_mps unset
// (Evolution::keep_mps keeps the evolved states). Measurements are taken at
// the measured steps as solve() takes them. The callback observes every
// iteration of the preparation and of the evolution; the evolution's
// iterations carry no energy (NaN).
ALPS_CLASSIC_DMRG_DECL EvolutionResult evolve(models::CompiledOperator const& initial,
    models::CompiledOperator const& final_hamiltonian, models::PhysicalBasis const& basis,
    Evolution const& evolution, Options const& options = {}, Measurements const& measurements = {},
    IterationCallback callback = {});

// The same for two models on one physical basis, each validated as solve()
// validates it.
ALPS_CLASSIC_DMRG_DECL EvolutionResult evolve(models::Hamiltonian const& initial,
    models::Hamiltonian const& final_hamiltonian, Evolution const& evolution, Options const& options = {},
    Measurements const& measurements = {}, IterationCallback callback = {});

// Measurement builders over a model's local operator names. `site_terms`
// holds one term per site that defines the operator (in site order);
// `pair_terms` holds one term per ordered site pair (i, j) over the sites
// that define both operators, i-major, with the on-site product at i == j.
// Both lower every pair separately, so no two columns are ever combined.
// A term that lowers to zero keeps its column as a scalar zero.
ALPS_CLASSIC_DMRG_DECL models::CompiledOperator site_terms(models::Hamiltonian const& model,
                                                           std::string const& op);
ALPS_CLASSIC_DMRG_DECL models::CompiledOperator pair_terms(models::Hamiltonian const& model,
                                                           std::string const& first,
                                                           std::string const& second);

} // namespace alps::dmrg::classic
