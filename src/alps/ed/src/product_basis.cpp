// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/assembly/product_basis.cpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names.
#include <alps/ed/product_basis.hpp>
#include "product_enumeration.hpp"
#include <stdexcept>
namespace alps::ed {
ProductBasis enumerate_product_basis(std::vector<int> const& site_dims,
                                     std::vector<BasisConstraint> const& constraints,
                                     std::size_t max_full_dimension) {
    ProductBasis result;
    struct Collector {
        ProductBasis &result;
        std::size_t sites;
        void reserve(std::size_t dimension) {
            if (sites && dimension > result.states.max_size() / sites)
                throw std::length_error("basis state buffer exceeds the native size limit");
            result.states.reserve(dimension * sites);
        }
        void operator()(std::vector<std::int32_t> const& state) {
            if (state.size() > result.states.max_size() - result.states.size())
                throw std::length_error("basis state buffer exceeds the native size limit");
            result.states.insert(result.states.end(), state.begin(), state.end());
        }
    };
    result.dimension = detail::visit_product_basis(site_dims, constraints, max_full_dimension,
                                                   Collector{result, site_dims.size()});
    return result;
}
} // namespace alps::ed
