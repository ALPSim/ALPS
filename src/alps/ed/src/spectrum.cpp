// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/runtime/spectrum.cpp.
// Adaptation: C++17 (no std::span, std::ranges or std::numbers), <alps/ed/...>
// include names, namespace alps::ed; the matrix and its blocks are held as
// canonical CSR (alpsx used Eigen sparse matrices) and solved by
// block_solver.hpp over the LAPACK and IETL drivers.
#include <alps/ed/spectrum.hpp>
#include <alps/ed/sparse_from_local.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include "block_solver.hpp"
#include "eigenstate_policy.hpp"
#include "exact_complex.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace alps::ed {
namespace {
using Complex = std::complex<double>;
using Matrix = SparseMatrix<Complex>;
constexpr double charge_tolerance = 1e-9;
double const pi = std::acos(-1.);

[[noreturn]] void invalid(const char* message) {
    throw std::invalid_argument(std::string("symmetry_spectrum: ") + message);
}

struct StateView {
    std::int64_t const* data;
    std::size_t size;
};
struct StateHash {
    std::size_t operator()(StateView state) const noexcept {
        std::size_t hash = 0;
        for (std::size_t i = 0; i < state.size; ++i)
            hash ^= std::hash<std::int64_t>{}(state.data[i]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        return hash;
    }
};
struct StateEqual {
    bool operator()(StateView a, StateView b) const noexcept {
        return a.size == b.size && std::equal(a.data, a.data + a.size, b.data);
    }
};
struct Action {
    std::size_t row;
    double sign;
};
struct Orbit {
    std::vector<Action> translated;
    std::vector<std::size_t> members;
    std::size_t charge;
};

void validate_matrix_input(std::vector<std::int64_t> const& indptr, std::vector<std::int64_t> const& indices,
                           std::vector<Complex> const& values, std::size_t dimension) {
    if (dimension > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || values.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        invalid("matrix exceeds the sparse-index range");
    if (indptr.size() != dimension + 1 || indices.size() != values.size())
        invalid("CSR array sizes do not match the matrix dimension");
    if (indptr.front() != 0 || indptr.back() != static_cast<std::int64_t>(values.size())
        || !std::is_sorted(indptr.begin(), indptr.end()))
        invalid("CSR offsets must increase from zero to the data length");
    for (std::size_t pos = 0; pos < values.size(); ++pos) {
        if (indices[pos] < 0 || indices[pos] >= static_cast<std::int64_t>(dimension))
            invalid("CSR column index is out of range");
        if (!std::isfinite(values[pos].real()) || !std::isfinite(values[pos].imag()))
            invalid("matrix values must be finite");
    }
}

Matrix read_matrix(std::vector<std::int64_t> const& indptr, std::vector<std::int64_t> const& indices,
                   std::vector<Complex> const& values, std::vector<std::size_t> const& row_charge) {
    auto const dimension = row_charge.size();
    CooFromLocal entries;
    entries.rows.reserve(values.size());
    entries.cols.reserve(values.size());
    entries.data.reserve(values.size());
    for (std::size_t row = 0; row < dimension; ++row) {
        std::map<std::int64_t, detail::ExactComplex> forbidden;
        for (auto pos = indptr[row]; pos < indptr[row + 1]; ++pos) {
            auto const p = static_cast<std::size_t>(pos);
            if (row_charge[row] != row_charge[static_cast<std::size_t>(indices[p])]) {
                // Check the original CSR entries, before duplicate reduction
                // can lose a residual or manufacture one. Allowed entries keep
                // their canonical accumulation and numerical tolerances.
                if (values[p] != Complex{}) forbidden[indices[p]] += detail::ExactComplex(values[p]);
            } else {
                entries.rows.push_back(static_cast<std::int64_t>(row));
                entries.cols.push_back(indices[p]);
                entries.data.push_back(values[p]);
            }
        }
        for (auto const& item : forbidden)
            if (!item.second.zero()) invalid("Hamiltonian does not conserve a declared charge");
        // Proven-zero forbidden entries never reach Hermiticity/translation
        // checks as a rounded duplicate-summation residue.
    }
    try {
        return canonical_csr<Complex>(entries, dimension, dimension);
    } catch (std::invalid_argument const&) {
        invalid("duplicate matrix-entry sums must be finite");
    }
}

Complex coefficient(Matrix const& matrix, std::size_t row, std::int64_t column) {
    auto const first = matrix.indices.begin() + matrix.offsets[row];
    auto const last = matrix.indices.begin() + matrix.offsets[row + 1];
    auto const it = std::lower_bound(first, last, column);
    if (it == last || *it != column) return Complex{};
    return matrix.values[static_cast<std::size_t>(it - matrix.indices.begin())];
}
} // namespace

std::vector<SpectrumBlock> symmetry_spectrum(
    std::vector<std::int64_t> const& indptr, std::vector<std::int64_t> const& indices,
    std::vector<Complex> const& values, std::size_t dimension, std::vector<std::int64_t> const& product_states,
    std::vector<SpectrumLocalSpace> const& local_spaces, std::vector<SpectrumTranslation> const& translations,
    int num_eigenvalues, std::vector<std::size_t> const& momentum_indices, std::string const& which,
    std::optional<double> tolerance, std::optional<int> max_iterations) {
    validate_eigensolver_controls(which, tolerance, max_iterations);
    if (num_eigenvalues <= 0) invalid("num_eigenvalues must be positive");
    std::size_t const nsites = local_spaces.size();
    if (nsites == 0 || product_states.size() / nsites != dimension || product_states.size() % nsites != 0)
        invalid("product states must have one column per local space");
    validate_matrix_input(indptr, indices, values, dimension);
    auto const nq = local_spaces.front().quantum_numbers.size();
    for (auto const& space : local_spaces) {
        if (space.parity.empty() || space.quantum_numbers.size() != nq)
            invalid("local spaces must have nonempty bases and the same charge labels");
        for (double parity : space.parity)
            if (parity != 1. && parity != -1.) invalid("local parity must be +1 or -1");
        for (auto const& charge : space.quantum_numbers) {
            if (charge.size() != space.parity.size()) invalid("local charge has the wrong dimension");
            for (double value : charge)
                if (!std::isfinite(value)) invalid("local charges must be finite");
        }
    }

    // Validate a representation of a product of cyclic translation groups.
    // This requires only primitive permutations, not an O(group_size^2)
    // multiplication table or a second geometry parser.
    std::vector<std::size_t> strides(translations.size());
    std::size_t group_size = 1;
    for (std::size_t d = translations.size(); d-- > 0;) {
        auto const& translation = translations[d];
        if (translation.period == 0 || translation.sites.size() != nsites)
            invalid("translation size or period is invalid");
        strides[d] = group_size;
        if (translation.period > std::numeric_limits<std::size_t>::max() / group_size)
            invalid("translation group is too large");
        group_size *= translation.period;
        std::vector<bool> seen(nsites);
        for (std::size_t site = 0; site < nsites; ++site) {
            auto const target = translation.sites[site];
            if (target >= nsites || seen[target]) invalid("translation must be a site permutation");
            seen[target] = true;
            if (local_spaces[site].parity != local_spaces[target].parity
                || local_spaces[site].quantum_numbers != local_spaces[target].quantum_numbers)
                invalid("translation must preserve each local space");
        }
        // Each cycle length must divide the declared period; do not loop
        // period times, since even a valid period can be arbitrarily large.
        std::fill(seen.begin(), seen.end(), false);
        std::size_t order = 1;
        for (std::size_t site = 0; site < nsites; ++site) {
            if (seen[site]) continue;
            std::size_t length = 0, target = site;
            do {
                seen[target] = true;
                target = translation.sites[target];
                ++length;
            } while (target != site);
            if (translation.period % length != 0) invalid("translation does not have the declared period");
            order = std::lcm(order, length); // Both factors divide period.
        }
        if (order != translation.period) invalid("translation period must equal its permutation order");
        for (std::size_t other = d + 1; other < translations.size(); ++other)
            for (std::size_t site = 0; site < nsites; ++site)
                if (translation.sites[translations[other].sites[site]]
                    != translations[other].sites[translation.sites[site]])
                    invalid("primitive translations must commute");
    }
    std::vector<std::size_t> momenta(momentum_indices);
    if (momenta.empty()) {
        momenta.resize(group_size);
        std::iota(momenta.begin(), momenta.end(), std::size_t{0});
    } else {
        std::sort(momenta.begin(), momenta.end());
        if (momenta.back() >= group_size || std::adjacent_find(momenta.begin(), momenta.end()) != momenta.end())
            invalid("requested momentum indices must be unique and within the translation grid");
    }

    std::unordered_map<StateView, std::size_t, StateHash, StateEqual> rows;
    rows.reserve(dimension);
    std::vector<std::vector<double>> totals(nq, std::vector<double>(dimension));
    for (std::size_t row = 0; row < dimension; ++row) {
        auto const* state = product_states.data() + row * nsites;
        for (std::size_t site = 0; site < nsites; ++site) {
            if (state[site] < 0 || static_cast<std::size_t>(state[site]) >= local_spaces[site].parity.size())
                invalid("product-state local index is out of range");
            for (std::size_t q = 0; q < nq; ++q)
                totals[q][row] += local_spaces[site].quantum_numbers[q][static_cast<std::size_t>(state[site])];
        }
        if (!rows.emplace(StateView{state, nsites}, row).second) invalid("product basis contains duplicate states");
    }
    // Cluster each charge against a fixed anchor. Approximate comparisons
    // must never serve as ordered-map comparators (they are nontransitive).
    std::vector<std::vector<std::size_t>> labels(dimension, std::vector<std::size_t>(nq));
    std::vector<std::vector<double>> charge_values(nq);
    for (std::size_t q = 0; q < nq; ++q) {
        std::vector<std::size_t> order(dimension);
        std::iota(order.begin(), order.end(), std::size_t{0});
        for (double total : totals[q])
            if (!std::isfinite(total)) invalid("total charge overflowed");
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return totals[q][a] < totals[q][b]; });
        for (auto row : order) {
            double const total = totals[q][row];
            if (charge_values[q].empty() || total - charge_values[q].back() > charge_tolerance)
                charge_values[q].push_back(total);
            labels[row][q] = charge_values[q].size() - 1;
        }
    }
    std::map<std::vector<std::size_t>, std::size_t> charge_ids;
    for (auto const& label : labels) charge_ids.try_emplace(label, 0);
    std::vector<std::vector<double>> charges;
    for (auto& item : charge_ids) {
        item.second = charges.size();
        std::vector<double> charge(nq);
        for (std::size_t q = 0; q < nq; ++q) charge[q] = charge_values[q][item.first[q]];
        charges.push_back(std::move(charge));
    }
    std::vector<std::size_t> row_charge(dimension);
    for (std::size_t row = 0; row < dimension; ++row) row_charge[row] = charge_ids.at(labels[row]);
    Matrix const matrix = read_matrix(indptr, indices, values, row_charge);

    std::vector<std::vector<Action>> generators(translations.size(), std::vector<Action>(dimension));
    std::vector<std::int64_t> translated(nsites);
    for (std::size_t d = 0; d < translations.size(); ++d) {
        auto const& permutation = translations[d].sites;
        for (std::size_t row = 0; row < dimension; ++row) {
            auto const* state = product_states.data() + row * nsites;
            double sign = 1.;
            for (std::size_t site = 0; site < nsites; ++site) {
                translated[permutation[site]] = state[site];
                if (local_spaces[site].parity[static_cast<std::size_t>(state[site])] > 0) continue;
                for (std::size_t other = site + 1; other < nsites; ++other)
                    if (permutation[site] > permutation[other]
                        && local_spaces[other].parity[static_cast<std::size_t>(state[other])] < 0)
                        sign = -sign;
            }
            auto const found = rows.find(StateView{translated.data(), nsites});
            if (found == rows.end()) invalid("product basis is not closed under translation");
            generators[d][row] = {found->second, sign};
            if (row_charge[row] != row_charge[found->second]) invalid("translation changes a conserved charge");
        }
    }

    double scale = 0.;
    for (auto z : matrix.values) scale = std::max({scale, std::abs(z.real()), std::abs(z.imag())});
    double const symmetry_tolerance = 1e-11 * scale;
    for (std::size_t row = 0; row < dimension; ++row)
        for (auto p = matrix.offsets[row]; p < matrix.offsets[row + 1]; ++p) {
            auto const column = matrix.indices[static_cast<std::size_t>(p)];
            auto const value = matrix.values[static_cast<std::size_t>(p)];
            if (std::abs(value - std::conj(coefficient(matrix, static_cast<std::size_t>(column), static_cast<std::int64_t>(row))))
                > symmetry_tolerance)
                invalid("Hamiltonian must be Hermitian");
            for (auto const& generator : generators) {
                auto const a = generator[row], b = generator[static_cast<std::size_t>(column)];
                if (std::abs(coefficient(matrix, a.row, static_cast<std::int64_t>(b.row)) - a.sign * b.sign * value)
                    > symmetry_tolerance)
                    invalid("Hamiltonian is not invariant under the lattice translations");
            }
        }
    if (dimension == 0) return {};

    std::vector<std::size_t> parent(group_size), step(group_size);
    for (std::size_t g = 1; g < group_size; ++g)
        for (std::size_t d = translations.size(); d-- > 0;)
            if ((g / strides[d]) % translations[d].period != 0) {
                parent[g] = g - strides[d];
                step[g] = d;
                break;
            }
    std::vector<Orbit> orbits;
    std::vector<bool> seen(dimension);
    for (std::size_t row = 0; row < dimension; ++row) {
        if (seen[row]) continue;
        Orbit orbit{{}, {}, row_charge[row]};
        orbit.translated.resize(group_size);
        orbit.translated[0] = {row, 1.};
        for (std::size_t g = 1; g < group_size; ++g) {
            auto const previous = orbit.translated[parent[g]];
            auto const action = generators[step[g]][previous.row];
            orbit.translated[g] = {action.row, previous.sign * action.sign};
        }
        for (auto action : orbit.translated) orbit.members.push_back(action.row);
        std::sort(orbit.members.begin(), orbit.members.end());
        orbit.members.erase(std::unique(orbit.members.begin(), orbit.members.end()), orbit.members.end());
        for (auto member : orbit.members) seen[member] = true;
        orbits.push_back(std::move(orbit));
    }

    std::vector<SpectrumBlock> result;
    std::size_t projected_dimension = 0;
    for (std::size_t k : momenta) {
        std::vector<double> momentum(translations.size());
        bool real_character = true;
        for (std::size_t d = 0; d < translations.size(); ++d) {
            auto const index = (k / strides[d]) % translations[d].period;
            momentum[d] = 2. * pi * static_cast<double>(index) / static_cast<double>(translations[d].period);
            real_character = real_character && (index == 0 || 2 * index == translations[d].period);
        }
        std::vector<Complex> character(group_size);
        for (std::size_t g = 0; g < group_size; ++g) {
            double phase = 0.;
            for (std::size_t d = 0; d < translations.size(); ++d)
                phase -= momentum[d] * static_cast<double>((g / strides[d]) % translations[d].period);
            character[g] = real_character ? Complex(std::cos(phase) < 0 ? -1. : 1.) : std::polar(1., phase);
        }
        std::vector<Complex> coefficients(dimension);
        std::vector<int> projected_row(dimension, -1);
        std::vector<std::size_t> sizes(charges.size());
        for (auto const& orbit : orbits) {
            for (std::size_t g = 0; g < group_size; ++g) {
                auto const action = orbit.translated[g];
                coefficients[action.row] += action.sign * character[g];
            }
            double norm = 0.;
            for (auto row : orbit.members) norm += std::norm(coefficients[row]);
            // A forbidden stabilizer character cancels exactly. Bound its
            // floating-point residue relative to the sum of group phases.
            double const roundoff = 64. * std::numeric_limits<double>::epsilon() * static_cast<double>(group_size);
            if (norm <= roundoff * roundoff * static_cast<double>(orbit.members.size())) continue;
            norm = std::sqrt(norm);
            for (auto row : orbit.members) {
                coefficients[row] /= norm;
                projected_row[row] = static_cast<int>(sizes[orbit.charge]);
            }
            ++sizes[orbit.charge];
        }
        std::vector<CooFromLocal> entries(charges.size());
        for (std::size_t row = 0; row < dimension; ++row) {
            if (projected_row[row] < 0) continue;
            auto& block = entries[row_charge[row]];
            for (auto p = matrix.offsets[row]; p < matrix.offsets[row + 1]; ++p) {
                auto const column = static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(p)]);
                if (projected_row[column] >= 0 && row_charge[row] == row_charge[column]) {
                    block.rows.push_back(projected_row[row]);
                    block.cols.push_back(projected_row[column]);
                    block.data.push_back(std::conj(coefficients[row]) * matrix.values[static_cast<std::size_t>(p)]
                                         * coefficients[column]);
                }
            }
        }
        for (std::size_t q = 0; q < charges.size(); ++q) {
            if (sizes[q] == 0) continue;
            auto const projected = canonical_csr<Complex>(entries[q], sizes[q], sizes[q]);
            auto solved = detail::solve_block(projected, which,
                detail::spectrum_plan(sizes[q], num_eigenvalues, which, tolerance, max_iterations), false);
            solved.spectrum.momentum = momentum;
            solved.spectrum.quantum_numbers = charges[q];
            result.push_back(std::move(solved.spectrum));
            projected_dimension += sizes[q];
        }
    }
    if (momenta.size() == group_size && projected_dimension != dimension)
        throw std::runtime_error("symmetry_spectrum: projected basis does not span the input basis");
    return result;
}
} // namespace alps::ed
