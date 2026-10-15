// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/runtime/finite_group.cpp.
// Adaptation: C++17 (no std::span or std::numbers), <alps/ed/...> include
// names, namespace alps::ed; the small representation matrices are handled
// by local row-major helpers (alpsx used Eigen).
#include <alps/ed/finite_group.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace alps::ed {
namespace {
using C = std::complex<double>;
using Dense = std::vector<C>; // row-major d x d
constexpr double tolerance = 1e-10;
double const pi = std::acos(-1.);
[[noreturn]] void invalid(char const* text) {
    throw std::invalid_argument(std::string("finite-group ED: ") + text);
}
Dense identity(std::size_t d) {
    Dense m(d * d);
    for (std::size_t i = 0; i < d; ++i) m[i * d + i] = 1.;
    return m;
}
Dense multiply(Dense const& a, Dense const& b, std::size_t d) {
    Dense m(d * d);
    for (std::size_t r = 0; r < d; ++r)
        for (std::size_t k = 0; k < d; ++k)
            for (std::size_t c = 0; c < d; ++c) m[r * d + c] += a[r * d + k] * b[k * d + c];
    return m;
}
Dense adjoint(Dense const& a, std::size_t d) {
    Dense m(d * d);
    for (std::size_t r = 0; r < d; ++r)
        for (std::size_t c = 0; c < d; ++c) m[r * d + c] = std::conj(a[c * d + r]);
    return m;
}
bool near(Dense const& a, Dense const& b) {
    for (std::size_t i = 0; i < a.size(); ++i)
        if (std::abs(a[i] - b[i]) > tolerance) return false;
    return true;
}
C trace(Dense const& a, std::size_t d) {
    C result{};
    for (std::size_t i = 0; i < d; ++i) result += a[i * d + i];
    return result;
}
} // namespace

struct FiniteGroupIrrep::Data {
    std::size_t dimension;
    std::string label;
    std::vector<std::vector<std::size_t>> generators;
    std::vector<std::vector<std::size_t>> permutations;
    std::vector<std::vector<C>> matrices;
    std::vector<C> characters;
};

FiniteGroupIrrep::FiniteGroupIrrep(std::vector<std::vector<std::size_t>> generators,
                                   std::vector<std::vector<C>> matrices, std::size_t dimension, std::string label,
                                   std::size_t max_order) {
    if (generators.empty() || generators.size() != matrices.size()) invalid("supply matching generators and matrices");
    if (!dimension || !max_order || dimension > max_order / dimension
        || dimension > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        invalid("invalid representation dimension or order limit");
    if (label.empty()) invalid("irrep label must not be empty");
    std::size_t const n = generators.front().size();
    Dense const one = identity(dimension);
    for (std::size_t g = 0; g < generators.size(); ++g) {
        if (generators[g].size() != n || matrices[g].size() != dimension * dimension)
            invalid("generator dimensions disagree");
        std::vector<bool> seen(n);
        for (auto site : generators[g]) {
            if (site >= n || seen[site]) invalid("site action must be a permutation");
            seen[site] = true;
        }
        for (auto z : matrices[g])
            if (!std::isfinite(z.real()) || !std::isfinite(z.imag())) invalid("representation entries must be finite");
        if (!near(multiply(adjoint(matrices[g], dimension), matrices[g], dimension), one))
            invalid("generator representation must be unitary");
    }
    auto data = std::make_shared<Data>();
    data->dimension = dimension;
    data->label = std::move(label);
    data->generators = generators;
    std::vector<std::size_t> id(n);
    std::iota(id.begin(), id.end(), std::size_t{0});
    data->permutations.push_back(id);
    data->matrices.push_back(one);
    std::map<std::vector<std::size_t>, std::size_t> indices{{id, 0}};
    for (std::size_t parent = 0; parent < data->permutations.size(); ++parent)
        for (std::size_t g = 0; g < generators.size(); ++g) {
            std::vector<std::size_t> permutation(n);
            for (std::size_t i = 0; i < n; ++i) permutation[i] = generators[g][data->permutations[parent][i]];
            Dense product = multiply(matrices[g], data->matrices[parent], dimension);
            auto found = indices.find(permutation);
            if (found != indices.end()) {
                if (!near(product, data->matrices[found->second])) invalid("representation violates a group relation");
            } else {
                if (data->permutations.size() == max_order) invalid("group closure exceeds max_order");
                indices.emplace(permutation, data->permutations.size());
                data->permutations.push_back(std::move(permutation));
                data->matrices.push_back(std::move(product));
            }
        }
    double norm = 0.;
    for (auto const& entry : data->matrices) {
        if (!near(multiply(adjoint(entry, dimension), entry, dimension), one))
            invalid("closed representation lost unitarity");
        auto const chi = trace(entry, dimension);
        data->characters.push_back(chi);
        norm += std::norm(chi);
    }
    if (std::abs(norm / static_cast<double>(data->permutations.size()) - 1.) > tolerance)
        invalid("representation is reducible");
    data_ = std::move(data);
}
std::size_t FiniteGroupIrrep::order() const { return data_->permutations.size(); }
std::size_t FiniteGroupIrrep::num_sites() const { return data_->permutations.front().size(); }
std::size_t FiniteGroupIrrep::dimension() const { return data_->dimension; }
std::string const& FiniteGroupIrrep::label() const { return data_->label; }
std::vector<std::vector<std::size_t>> const& FiniteGroupIrrep::generators() const { return data_->generators; }
std::vector<std::vector<std::size_t>> const& FiniteGroupIrrep::permutations() const { return data_->permutations; }
std::vector<std::vector<C>> const& FiniteGroupIrrep::matrices() const { return data_->matrices; }
std::vector<C> const& FiniteGroupIrrep::characters() const { return data_->characters; }

FiniteGroupIrrep dihedral_irrep(std::size_t n, std::size_t k, int reflection) {
    if (n < 3 || n > 2048 || k > n / 2 || (reflection != 1 && reflection != -1))
        invalid("invalid dihedral size, angular index or reflection sign");
    std::vector<std::vector<std::size_t>> generators(2, std::vector<std::size_t>(n));
    for (std::size_t i = 0; i < n; ++i) {
        generators[0][i] = (i + 1) % n;
        generators[1][i] = (n - i) % n;
    }
    bool const scalar = k == 0 || (n % 2 == 0 && k == n / 2);
    if (scalar)
        return FiniteGroupIrrep(generators, {{k == 0 ? 1. : -1.}, {static_cast<double>(reflection)}}, 1,
                                std::string(k == 0 ? "A" : "B") + (reflection == 1 ? "1" : "2"));
    if (reflection != 1) invalid("reflection sign selects only one-dimensional dihedral irreps");
    double const angle = 2. * pi * static_cast<double>(k) / static_cast<double>(n), c = std::cos(angle), s = std::sin(angle);
    return FiniteGroupIrrep(generators, {{c, -s, s, c}, {1., 0., 0., -1.}}, 2, "E" + std::to_string(k));
}

void validate_irrep_family(std::vector<FiniteGroupIrrep> const& irreps) {
    if (irreps.empty()) invalid("spectrum requires a complete irrep family");
    auto const& first = irreps.front();
    std::map<std::vector<std::size_t>, std::size_t> indices;
    for (std::size_t g = 0; g < first.order(); ++g) indices.emplace(first.permutations()[g], g);
    std::vector<std::vector<C>> characters;
    std::set<std::string> labels;
    std::size_t dimensions = 0;
    for (auto const& irrep : irreps) {
        if (irrep.order() != first.order() || !labels.insert(irrep.label()).second)
            invalid("irrep family has different groups or duplicate labels");
        if (irrep.dimension() > first.order() / irrep.dimension()) invalid("invalid irrep dimension");
        auto const square = irrep.dimension() * irrep.dimension();
        if (square > first.order() - dimensions) invalid("irrep family is overcomplete");
        dimensions += square;
        auto& row = characters.emplace_back(first.order());
        for (std::size_t g = 0; g < irrep.order(); ++g) {
            auto found = indices.find(irrep.permutations()[g]);
            if (found == indices.end()) invalid("irrep family must describe the same permutation group");
            row[found->second] = irrep.characters()[g];
        }
        for (std::size_t previous = 0; previous + 1 < characters.size(); ++previous) {
            C inner{};
            for (std::size_t g = 0; g < first.order(); ++g) inner += std::conj(characters[previous][g]) * row[g];
            if (std::abs(inner) / static_cast<double>(first.order()) > tolerance)
                invalid("irrep family contains equivalent representations");
        }
    }
    if (dimensions != first.order()) invalid("irrep family is incomplete");
}

} // namespace alps::ed
