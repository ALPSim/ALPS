// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/runtime/eigenstate_policy.hpp.
// Adaptation: <alps/ed/...> include names; the Lanczos threshold and defaults
// follow the IETL-based driver (8*count + 200 Krylov room, driver defaults).
#pragma once

#include <alps/ed/eigensolver.hpp>
#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace alps::ed::detail {

// Resolved dense/Lanczos choice and Lanczos controls for one block. Callers
// validate the public controls first; dense solves ignore the controls.
// Zero controls select the lanczos_csr defaults.
struct SolvePlan {
    std::size_t count = 0; // Requested count limited by the block dimension.
    bool iterative = false;
    double tolerance = 0.;
    int max_iterations = 0;
};

inline void validate_request(EigenstateOptions const& options) {
    validate_eigensolver_controls(options.which, options.tolerance, options.max_iterations);
    if (options.count <= 0) throw std::invalid_argument("eigenvalue count must be positive");
}

inline bool lanczos_requested(EigenstateMethod method, std::size_t dimension) {
    switch (method) {
    case EigenstateMethod::automatic: return dimension > 1000;
    case EigenstateMethod::dense: return false;
    case EigenstateMethod::lanczos: return true;
    }
    throw std::invalid_argument("unknown eigenstate method");
}

inline SolvePlan eigenstate_plan(EigenstateOptions const& options, std::size_t dimension) {
    SolvePlan plan;
    plan.count = std::min(dimension, static_cast<std::size_t>(options.count));
    plan.iterative = lanczos_requested(options.method, dimension) && dimension > 8 * plan.count + 200
        && !(options.method == EigenstateMethod::automatic && options.which == "SM");
    plan.tolerance = options.tolerance.value_or(0.);
    plan.max_iterations = options.max_iterations.value_or(0);
    return plan;
}

} // namespace alps::ed::detail
