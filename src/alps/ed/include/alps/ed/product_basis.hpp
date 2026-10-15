// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/product_basis.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL exports.

#pragma once

#include <alps/ed_export.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace alps::ed {

struct BasisConstraint {
    // values[site][local state]; every site must supply this quantum number.
    std::vector<std::vector<double>> values;
    double target = 0.;
};

struct ProductBasis {
    // Row-major (dimension, number of sites), in lexicographic order.
    // Site 0 is the most significant digit. The vacuum is one empty row.
    std::vector<std::int32_t> states;
    std::size_t dimension = 0;
};

// Enumerate only the selected sector, using bounds on remaining charges to
// prune branches. Charge sums accumulate from site 0 to the last site; every
// final charge must differ from its target by strictly less than 1e-9.
// No full-space integer encoding or recursive call stack is needed.
// max_full_dimension applies only when constraints is empty.
ALPS_ED_DECL ProductBasis enumerate_product_basis(
    std::vector<int> const& site_dims,
    std::vector<BasisConstraint> const& constraints = {},
    std::size_t max_full_dimension = 10'000'000);

} // namespace alps::ed
