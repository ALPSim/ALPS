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
