// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The algebraic QMC sign-problem check of the legacy <alps/model/sign.h>,
// on compiled operators (implementation and provenance in src/sign.cpp).
#pragma once

#include <alps/models_export.h>
#include <alps/models/compiled.hpp>
#include <alps/models/hamiltonian.hpp>

namespace alps {
namespace models {

// Whether the sublattice-sign rule leaves positive off-diagonal weights of
// -H: the one-site terms are summed per site and the two-site terms per
// site pair, each sum must have off-diagonal elements of one sign (bond
// elements count where both sites change), and the bond and site signs must
// admit a sublattice assignment (the legacy frustration test). Elements
// within 1e-10 of zero are absent, as in the legacy check. Unlike it, every
// interacting pair is checked on its own, so inhomogeneous couplings and
// pairs off the lattice take part. Purely algebraic: no sampling.
//
// Throws std::invalid_argument for operators outside the rule: terms on more
// than two sites, Jordan-Wigner strings (fermion hopping), and matrix
// elements with an imaginary part above 1e-10.
ALPS_MODELS_DECL bool has_sign_problem(CompiledOperator const& op);
// has_sign_problem(compile_hamiltonian(model).op).
ALPS_MODELS_DECL bool has_sign_problem(Hamiltonian const& model);

} // namespace models
} // namespace alps
