// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The validation compiled-operator entry points share (defined in
// projected_problem.cpp).
#pragma once

#include <alps/ed/projected_problem.hpp>

namespace alps::ed::detail {

// Hermiticity of the full operator and conservation of every charge of
// `spaces` (sorted, unique names; one row per local space), with the operator's
// Jordan-Wigner strings materialized as diagonal factors. Throws
// invalid_argument on a violation or a malformed operator.
void validate_compiled_model(models::CompiledOperator const& op, ProjectedSpaces const& spaces);

} // namespace alps::ed::detail
