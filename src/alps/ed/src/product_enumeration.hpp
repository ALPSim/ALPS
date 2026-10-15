// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/assembly/product_enumeration.hpp.
// Adaptation: C++17 (no std::span or requires-expressions), <alps/ed/...> include names.
#pragma once
#include <alps/ed/product_basis.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
namespace alps::ed::detail {

template <class Visitor, class = void> struct has_reserve : std::false_type {};
template <class Visitor>
struct has_reserve<Visitor, std::void_t<decltype(std::declval<Visitor&>().reserve(std::size_t{}))>>
    : std::true_type {};

// Shared traversal; the visitor borrows one scratch row only until it returns.
template <class Visitor>
std::size_t visit_product_basis(std::vector<int> const& site_dims, std::vector<BasisConstraint> const& constraints,
                                std::size_t max_full_dimension, Visitor &&visit) {
    constexpr double tolerance = 1e-9;
    auto const sites = site_dims.size(), charges = constraints.size();
    for (int dim : site_dims)
        if (dim < 1)
            throw std::invalid_argument("basis site dimensions must be positive");
    for (auto const &constraint : constraints) {
        if (!std::isfinite(constraint.target))
            throw std::invalid_argument("basis sector targets must be finite");
        if (constraint.values.size() != sites)
            throw std::invalid_argument("basis charge site count mismatch");
        for (std::size_t site = 0; site < sites; ++site) {
            auto const &values = constraint.values[site];
            if (values.size() != static_cast<std::size_t>(site_dims[site]))
                throw std::invalid_argument("basis charge values must match local dimensions");
            if (!std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); }))
                throw std::invalid_argument("basis charge values must be finite");
        }
    }
    std::size_t count = 0;
    std::vector<std::int32_t> state(sites, 0);
    if (constraints.empty()) {
        std::size_t dimension = 1;
        for (int dim : site_dims) {
            if (dimension > max_full_dimension / static_cast<std::size_t>(dim))
                throw std::invalid_argument(
                    "Basis dimension too large for full enumeration; use a QN sector");
            dimension *= static_cast<std::size_t>(dim);
        }
        if (dimension > max_full_dimension)
            throw std::invalid_argument("Basis dimension too large for full enumeration; use a QN sector");
        if constexpr (has_reserve<std::decay_t<Visitor>>::value)
            visit.reserve(dimension);

        for (std::size_t row = 0; row < dimension; ++row) {
            visit(state);
            for (std::size_t site = sites; site-- > 0;) {
                if (++state[site] < site_dims[site])
                    break;
                state[site] = 0;
            }
        }
        return dimension;
    }

    // Suffix bounds use the same right-to-left accumulation as the Python
    // reference. Partial charges below still accumulate left-to-right.
    std::vector<std::vector<double>> remaining_min(sites + 1, std::vector<double>(charges, 0.));
    auto remaining_max = remaining_min;
    for (std::size_t site = sites; site-- > 0;)
        for (std::size_t q = 0; q < charges; ++q) {
            auto const &values = constraints[q].values[site];
            auto const [low, high] = std::minmax_element(values.begin(), values.end());
            remaining_min[site][q] = remaining_min[site + 1][q] + *low;
            remaining_max[site][q] = remaining_max[site + 1][q] + *high;
        }
    std::vector<std::vector<double>> partial(sites + 1, std::vector<double>(charges, 0.));
    std::vector<int> next_local(sites, 0);
    std::size_t depth = 0;
    while (true) {
        if (depth == sites) {
            bool matches = true;
            for (std::size_t q = 0; q < charges; ++q)
                matches = matches && std::abs(partial[depth][q] - constraints[q].target) < tolerance;
            if (matches) {
                if (count == std::numeric_limits<std::size_t>::max())
                    throw std::length_error("basis dimension exceeds the native size limit");
                visit(state);
                ++count;
            }
            if (depth == 0)
                break;
            --depth;
            continue;
        }
        if (next_local[depth] == site_dims[depth]) {
            next_local[depth] = 0;
            if (depth == 0)
                break;
            --depth;
            continue;
        }
        int const local = next_local[depth]++;
        state[depth] = local;
        bool reachable = true;
        for (std::size_t q = 0; q < charges; ++q) {
            double const sum = partial[depth][q] + constraints[q].values[depth][static_cast<std::size_t>(local)];
            partial[depth + 1][q] = sum;
            if (depth + 1 < sites) {
                double const needed = constraints[q].target - sum;
                if (needed < remaining_min[depth + 1][q] - tolerance ||
                    needed > remaining_max[depth + 1][q] + tolerance)
                    reachable = false;
            } else if (std::abs(sum - constraints[q].target) > tolerance) {
                reachable = false;
            }
        }
        if (reachable)
            ++depth;
    }
    return count;
}
} // namespace alps::ed::detail
