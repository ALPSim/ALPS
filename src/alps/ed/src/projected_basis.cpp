// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/runtime/projected_basis.cpp.
// Adaptation: C++17 (no std::span, templated lambdas, std::erase_if,
// map::contains, std::ranges or std::numbers), <alps/ed/...> include names,
// namespace alps::ed; blocks are held as canonical CSR (alpsx used Eigen
// sparse matrices) and solved by block_solver.hpp over the LAPACK and IETL
// drivers; the finite point-group basis is not ported.
#include <alps/ed/projected_basis.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/ed/eigensolver.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include "block_solver.hpp"
#include "eigenstate_policy.hpp"
#include "exact_complex.hpp"
#include "local_action.hpp"
#include "product_enumeration.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace alps::ed {
namespace detail {
using Complex = std::complex<double>;
using State = std::vector<std::int32_t>;
struct ProjectedData {
    std::vector<SpectrumLocalSpace> spaces;
    std::vector<SpectrumTranslation> translations;
    std::vector<int> dims;
    std::vector<std::optional<double>> targets;
    std::vector<BasisConstraint> constraints;
    std::vector<double> momentum;
    std::size_t momentum_index = 0, parent_dimension = 0;
    std::vector<std::vector<std::size_t>> permutations;
    std::vector<Complex> characters;
    std::vector<std::int32_t> representatives;
    std::vector<std::size_t> sizes;
    std::map<State, std::size_t> indices;
    std::vector<std::vector<double>> charge_anchors;
    bool graded = false;
};
struct ProjectedAccess {
    static ProjectedData const& get(ProjectedBasis const& basis) { return *basis.data_; }
    static ProjectedBasis make(std::shared_ptr<const ProjectedData> data) {
        return ProjectedBasis(std::move(data));
    }
};
} // namespace detail
namespace {
using detail::Complex;
using detail::ProjectedAccess;
using detail::ProjectedData;
using detail::State;
using Action = std::map<State, Complex>;
using Matrix = SparseMatrix<Complex>;
using Prepared = detail::PreparedTerm<Complex>;
double const pi = std::acos(-1.);
[[noreturn]] void invalid(char const* text) {
    throw std::invalid_argument(std::string("projected ED: ") + text);
}
bool finite(Complex z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }

std::vector<double> charges(ProjectedData const& d, State const& state) {
    std::vector<double> result(d.targets.size());
    for (std::size_t site = 0; site < state.size(); ++site)
        for (std::size_t q = 0; q < result.size(); ++q)
            result[q] += d.spaces[site].quantum_numbers[q][static_cast<std::size_t>(state[site])];
    for (double q : result)
        if (!std::isfinite(q)) invalid("total charge overflowed");
    return result;
}
bool selected(ProjectedData const& d, State const& state) {
    auto totals = charges(d, state);
    for (std::size_t q = 0; q < totals.size(); ++q)
        if (d.targets[q] && std::abs(totals[q] - *d.targets[q]) >= 1e-9) return false;
    return true;
}
std::vector<std::size_t> charge_labels(ProjectedData const& d, std::vector<double> const& totals) {
    std::vector<std::size_t> result;
    for (std::size_t q = 0; q < totals.size(); ++q) {
        auto const& anchors = d.charge_anchors[q];
        auto it = std::upper_bound(anchors.begin(), anchors.end(), totals[q]);
        if (it == anchors.begin() || totals[q] - *(it - 1) > 1e-9)
            result.push_back(std::numeric_limits<std::size_t>::max());
        else
            result.push_back(static_cast<std::size_t>(it - anchors.begin() - 1));
    }
    return result;
}
double translate(ProjectedData const& d, State const& state, std::vector<std::size_t> const& permutation,
                 State& out) {
    out.resize(state.size());
    double sign = 1.;
    for (std::size_t site = 0; site < state.size(); ++site) {
        out[permutation[site]] = state[site];
        if (!d.graded || d.spaces[site].parity[static_cast<std::size_t>(state[site])] > 0) continue;
        for (std::size_t other = site + 1; other < state.size(); ++other)
            if (permutation[site] > permutation[other]
                && d.spaces[other].parity[static_cast<std::size_t>(state[other])] < 0)
                sign = -sign;
    }
    return sign;
}
struct Canonical {
    State state;
    std::size_t group;
    double sign;
};
Canonical canonical(ProjectedData const& d, State const& state) {
    State smallest(state), scratch;
    std::size_t group = 0;
    double phase_sign = 1.;
    for (std::size_t g = 1; g < d.permutations.size(); ++g) {
        double const sign = translate(d, state, d.permutations[g], scratch);
        if (scratch < smallest) {
            smallest = scratch;
            group = g;
            phase_sign = sign;
        }
    }
    return {std::move(smallest), group, phase_sign};
}
std::optional<std::pair<std::size_t, Complex>> project(ProjectedData const& d, State const& state) {
    auto c = canonical(d, state);
    auto it = d.indices.find(c.state);
    if (it == d.indices.end()) return {};
    return std::pair<std::size_t, Complex>{it->second, c.sign * std::conj(d.characters[c.group])
                                                           / std::sqrt(static_cast<double>(d.sizes[it->second]))};
}
State row(ProjectedData const& d, std::size_t index) {
    auto const first = d.representatives.begin() + static_cast<std::ptrdiff_t>(index * d.dims.size());
    return State(first, first + static_cast<std::ptrdiff_t>(d.dims.size()));
}
std::map<State, Complex> orbit(ProjectedData const& d, std::size_t index) {
    auto const rep = row(d, index);
    std::map<State, Complex> result;
    State scratch;
    for (std::size_t g = 0; g < d.permutations.size(); ++g) {
        double const sign = translate(d, rep, d.permutations[g], scratch);
        result.try_emplace(scratch, sign * d.characters[g] / std::sqrt(static_cast<double>(d.sizes[index])));
    }
    return result;
}
template<class Emit>
void project_components(ProjectedData const& d, State const& state, Emit&& emit) {
    if (auto p = project(d, state)) emit(p->first, p->second);
}

std::shared_ptr<ProjectedData> metadata(std::vector<SpectrumLocalSpace> const& spaces,
                                        std::vector<SpectrumTranslation> const& translations,
                                        std::vector<std::optional<double>> const& targets,
                                        std::size_t momentum_index) {
    auto d = std::make_shared<ProjectedData>();
    d->spaces = spaces;
    d->translations = translations;
    auto const nq = spaces.empty() ? 0 : spaces.front().quantum_numbers.size(), n = spaces.size();
    if (!targets.empty() && targets.size() != nq) invalid("charge target count mismatch");
    d->targets.resize(nq);
    if (!targets.empty()) std::copy(targets.begin(), targets.end(), d->targets.begin());
    for (auto target : d->targets)
        if (target && !std::isfinite(*target)) invalid("charge targets must be finite");
    for (auto const& space : spaces) {
        if (space.parity.empty()
            || space.parity.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
            || space.quantum_numbers.size() != nq)
            invalid("invalid local space dimensions");
        d->dims.push_back(static_cast<int>(space.parity.size()));
        for (double p : space.parity) {
            if (p != 1. && p != -1.) invalid("local parity must be +1 or -1");
            d->graded = d->graded || p < 0;
        }
        for (auto const& q : space.quantum_numbers) {
            if (q.size() != space.parity.size()) invalid("local charge dimension mismatch");
            for (double x : q)
                if (!std::isfinite(x)) invalid("local charges must be finite");
        }
    }
    for (std::size_t q = 0; q < nq; ++q)
        if (d->targets[q]) {
            BasisConstraint c;
            c.target = *d->targets[q];
            for (auto const& space : spaces) c.values.push_back(space.quantum_numbers[q]);
            d->constraints.push_back(std::move(c));
        }
    std::size_t group_size = 1;
    std::vector<std::size_t> strides(translations.size());
    for (std::size_t direction = translations.size(); direction-- > 0;) {
        auto const& t = translations[direction];
        if (!t.period || t.sites.size() != n) invalid("translation size or period is invalid");
        strides[direction] = group_size;
        if (t.period > std::numeric_limits<std::size_t>::max() / group_size)
            invalid("translation group is too large");
        group_size *= t.period;
        std::vector<bool> seen(n);
        for (std::size_t site = 0; site < n; ++site) {
            auto const dst = t.sites[site];
            if (dst >= n || seen[dst]) invalid("translation must be a permutation");
            seen[dst] = true;
            if (spaces[site].parity != spaces[dst].parity
                || spaces[site].quantum_numbers != spaces[dst].quantum_numbers)
                invalid("translation must preserve local spaces");
        }
        std::fill(seen.begin(), seen.end(), false);
        std::size_t order = 1;
        for (std::size_t site = 0; site < n; ++site)
            if (!seen[site]) {
                std::size_t length = 0, at = site;
                do {
                    seen[at] = true;
                    at = t.sites[at];
                    ++length;
                } while (at != site);
                if (t.period % length) invalid("translation does not have the declared period");
                order = std::lcm(order, length);
            }
        if (order != t.period) invalid("translation period must equal its permutation order");
        for (std::size_t other = direction + 1; other < translations.size(); ++other)
            for (std::size_t site = 0; site < n; ++site)
                if (t.sites[translations[other].sites[site]] != translations[other].sites[t.sites[site]])
                    invalid("primitive translations must commute");
    }
    if (momentum_index >= group_size) invalid("momentum index outside translation grid");
    d->momentum_index = momentum_index;
    d->momentum.resize(translations.size());
    bool real = true;
    for (std::size_t direction = 0; direction < translations.size(); ++direction) {
        auto const k = (momentum_index / strides[direction]) % translations[direction].period;
        d->momentum[direction] = 2. * pi * static_cast<double>(k) / static_cast<double>(translations[direction].period);
        real = real && (k == 0 || (translations[direction].period % 2 == 0 && k == translations[direction].period / 2));
    }
    d->permutations.resize(group_size, std::vector<std::size_t>(n));
    d->characters.resize(group_size);
    std::iota(d->permutations[0].begin(), d->permutations[0].end(), std::size_t{0});
    for (std::size_t g = 0; g < group_size; ++g) {
        double phase = 0.;
        for (std::size_t direction = 0; direction < translations.size(); ++direction)
            phase -= d->momentum[direction]
                * static_cast<double>((g / strides[direction]) % translations[direction].period);
        d->characters[g] = real ? Complex(std::cos(phase) < 0 ? -1. : 1.) : std::polar(1., phase);
        if (g)
            for (std::size_t direction = translations.size(); direction-- > 0;)
                if ((g / strides[direction]) % translations[direction].period) {
                    auto const parent = g - strides[direction];
                    for (std::size_t site = 0; site < n; ++site)
                        d->permutations[g][site] = translations[direction].sites[d->permutations[parent][site]];
                    break;
                }
    }
    return d;
}

void fill_basis(ProjectedData& d, bool all_orbits = false) {
    std::vector<std::set<double>> unique(d.targets.size());
    State scratch;
    d.parent_dimension = detail::visit_product_basis(
        d.dims, d.constraints, std::numeric_limits<std::size_t>::max(), [&](State const& state) {
            auto totals = charges(d, state);
            for (std::size_t q = 0; q < totals.size(); ++q) unique[q].insert(totals[q]);
            for (std::size_t g = 1; g < d.permutations.size(); ++g) {
                translate(d, state, d.permutations[g], scratch);
                if (std::lexicographical_compare(scratch.begin(), scratch.end(), state.begin(), state.end())) {
                    if (!selected(d, scratch)) invalid("selected sector is not closed under translation");
                    return;
                }
            }
            Complex stabilizer_sum = 0.;
            std::size_t stabilizers = 0;
            for (std::size_t g = 0; g < d.permutations.size(); ++g) {
                double const sign = translate(d, state, d.permutations[g], scratch);
                if (!selected(d, scratch)) invalid("selected sector is not closed under translation");
                if (std::equal(scratch.begin(), scratch.end(), state.begin())) {
                    ++stabilizers;
                    stabilizer_sum += sign * d.characters[g];
                }
            }
            if (!all_orbits && std::abs(stabilizer_sum)
                    <= 64. * std::numeric_limits<double>::epsilon() * static_cast<double>(d.permutations.size()))
                return;
            auto const index = d.sizes.size();
            d.indices.emplace(state, index);
            d.representatives.insert(d.representatives.end(), state.begin(), state.end());
            d.sizes.push_back(d.permutations.size() / stabilizers);
        });
    d.charge_anchors.resize(unique.size());
    for (std::size_t q = 0; q < unique.size(); ++q)
        for (double value : unique[q])
            if (d.charge_anchors[q].empty() || value - d.charge_anchors[q].back() > 1e-9)
                d.charge_anchors[q].push_back(value);
}

struct SingleState {
    using State = detail::State;
    State const& input;
    int digit(std::size_t, int site) const { return input[static_cast<std::size_t>(site)]; }
    State state(std::size_t) const { return input; }
    void shift(State& target, int site, int, int after) const { target[static_cast<std::size_t>(site)] = after; }
};
struct Operator {
    std::vector<Prepared> terms, adjoint;
    models::CompiledOperator tensors;
    std::vector<std::vector<double> const*> no_parity;
};
Operator prepare_operator(models::CompiledOperator const& input, ProjectedData const& d) {
    if (input.site_dims != d.dims) invalid("operator local dimensions do not match basis");
    auto const checked = detail::prepare_terms<Complex>(
        input,
        [](Complex z, char const*) {
            if (!finite(z)) invalid("operator coefficients and entries must be finite");
            return z;
        },
        true);
    (void)checked;
    Operator result;
    result.tensors.site_dims = input.site_dims;
    result.no_parity.resize(d.dims.size(), nullptr);
    // Fold the input-state JW signs into tensor columns once. This also makes
    // Hermitian conjugation exact without commuting strings past local ops.
    for (auto const& term : input.terms) {
        std::map<int, std::vector<Complex>> factors;
        for (auto const& factor : term.factors) factors.emplace(factor.site, factor.matrix);
        for (int site : term.jordan_wigner) {
            auto const dim = static_cast<std::size_t>(d.dims[static_cast<std::size_t>(site)]);
            auto [it, created] = factors.try_emplace(site, dim * dim, Complex{});
            if (created)
                for (std::size_t i = 0; i < dim; ++i) it->second[i * dim + i] = 1.;
            for (std::size_t r = 0; r < dim; ++r)
                for (std::size_t c = 0; c < dim; ++c)
                    it->second[r * dim + c] *= d.spaces[static_cast<std::size_t>(site)].parity[c];
        }
        models::LocalTerm folded{term.coefficient, {}, {}};
        for (auto& [site, matrix] : factors) folded.factors.push_back({site, std::move(matrix)});
        result.tensors.terms.push_back(std::move(folded));
    }
    auto dagger = result.tensors;
    for (auto& term : dagger.terms) {
        term.coefficient = std::conj(term.coefficient);
        for (auto& factor : term.factors) {
            auto const dim = static_cast<std::size_t>(d.dims[static_cast<std::size_t>(factor.site)]);
            auto const old = factor.matrix;
            for (std::size_t r = 0; r < dim; ++r)
                for (std::size_t c = 0; c < dim; ++c) factor.matrix[r * dim + c] = std::conj(old[c * dim + r]);
        }
    }
    auto const convert = [](Complex z, char const*) { return z; };
    result.terms = detail::prepare_terms<Complex>(result.tensors, convert, true);
    result.adjoint = detail::prepare_terms<Complex>(dagger, convert, true);
    return result;
}
template<class Emit>
void act(Operator const& op, State const& state, Emit&& emit, bool adjoint = false) {
    SingleState basis{state};
    for (auto const& term : adjoint ? op.adjoint : op.terms) {
        if (term.coeff == Complex{}) continue;
        detail::apply_term_state(term, basis, op.no_parity, 0, [&](State const& target, Complex value) {
            if (!finite(value)) invalid("operator action overflowed");
            if (value != Complex{}) emit(target, value);
        });
    }
}
detail::ExactComplex exact_amplitude(Operator const& op, std::size_t term, State const& source,
                                     State const& target, bool adjoint) {
    auto const& input = op.tensors; // JW signs have already been folded exactly into columns.
    detail::ExactComplex value(adjoint ? std::conj(input.terms[term].coefficient) : input.terms[term].coefficient);
    for (auto const& factor : input.terms[term].factors) {
        auto const site = static_cast<std::size_t>(factor.site);
        auto const dim = static_cast<std::size_t>(input.site_dims[site]);
        auto const& matrix = factor.matrix;
        auto const s = static_cast<std::size_t>(source[site]), t = static_cast<std::size_t>(target[site]);
        auto const entry = adjoint ? std::conj(matrix[s * dim + t]) : matrix[t * dim + s];
        value *= detail::ExactComplex(entry);
    }
    return value;
}
Action action(Operator const& op, ProjectedData const& d, State const& state, bool adjoint = false,
              bool require_conservation = false) {
    struct Amplitude {
        Complex rounded{};
        std::optional<detail::ExactComplex> forbidden;
    };
    std::map<State, Amplitude> amplitudes;
    auto const source_labels = charge_labels(d, charges(d, state));
    SingleState basis{state};
    auto const& terms = adjoint ? op.adjoint : op.terms;
    for (std::size_t t = 0; t < terms.size(); ++t) {
        auto const& term = terms[t];
        if (term.coeff == Complex{}) continue;
        // The scatter enumerates nonzero local support even if an intermediate
        // double product underflows. Do not discard a rounded zero before the
        // forbidden transition has been reconstructed from its original factors.
        detail::apply_term_state(term, basis, op.no_parity, 0, [&](State const& target, Complex value) {
            auto [it, inserted] = amplitudes.try_emplace(target);
            auto& amplitude = it->second;
            if (inserted && (!selected(d, target) || charge_labels(d, charges(d, target)) != source_labels))
                amplitude.forbidden.emplace();
            if (amplitude.forbidden) *amplitude.forbidden += exact_amplitude(op, t, state, target, adjoint);
            else if (!finite(value)) invalid("operator action overflowed");
            if (value != Complex{}) amplitude.rounded += value;
        });
    }
    Action result;
    for (auto const& [target, amplitude] : amplitudes) {
        if (amplitude.forbidden) {
            // Exact cancellation must not leave a rounded residue that then
            // fabricates a Hermiticity or translation failure. Nonzero escapes
            // from the selected parent sector have no amplitude tolerance.
            if (amplitude.forbidden->zero()) continue;
            if (require_conservation) invalid("Hamiltonian does not conserve a declared charge");
        }
        if (!finite(amplitude.rounded)) invalid("summed operator action overflowed");
        result.emplace(target, amplitude.rounded);
    }
    return result;
}
using TensorKey = std::vector<std::pair<int, std::vector<double>>>;
std::map<TensorKey, Complex> tensor_sum(Operator const& op, ProjectedData const& d,
                                        std::vector<std::size_t> const* permutation = nullptr,
                                        bool dagger = false) {
    std::map<TensorKey, Complex> out;
    for (auto const& term : op.tensors.terms) {
        TensorKey key;
        for (auto const& factor : term.factors) {
            int const site = factor.site, dim = d.dims[static_cast<std::size_t>(site)];
            auto const& matrix = factor.matrix;
            bool identity = true;
            std::vector<double> values;
            for (int r = 0; r < dim; ++r)
                for (int c = 0; c < dim; ++c) {
                    auto const z = dagger ? std::conj(matrix[static_cast<std::size_t>(c * dim + r)])
                                          : matrix[static_cast<std::size_t>(r * dim + c)];
                    identity = identity && (z == Complex(r == c ? 1. : 0.));
                    values.push_back(z.real());
                    values.push_back(z.imag());
                }
            if (!identity)
                key.emplace_back(permutation ? static_cast<int>((*permutation)[static_cast<std::size_t>(site)]) : site,
                                 std::move(values));
        }
        std::sort(key.begin(), key.end());
        out[std::move(key)] += dagger ? std::conj(term.coefficient) : term.coefficient;
    }
    for (auto it = out.begin(); it != out.end();)
        it = it->second == Complex{} ? out.erase(it) : std::next(it);
    return out;
}
bool certified(Operator const& op, ProjectedData const& d) {
    if (d.graded) return false;
    auto const original = tensor_sum(op, d);
    if (original != tensor_sum(op, d, nullptr, true)) return false;
    for (auto const& t : d.translations)
        if (original != tensor_sum(op, d, &t.sites)) return false;
    // Fixed exact charge deltas prove conservation termwise. Unusual floating
    // charges and mixed-delta matrices use the general action check instead.
    for (std::size_t q = 0; q < d.targets.size(); ++q) {
        double bound = 0.;
        for (auto const& space : d.spaces) {
            double largest = 0.;
            for (double x : space.quantum_numbers[q]) {
                if (std::abs(x) > 1e12 || 2. * x != std::round(2. * x)) return false;
                largest = std::max(largest, std::abs(x));
            }
            bound += largest;
        }
        if (bound > 0x1p51) return false; // Every half-integer partial sum is exact.
    }
    for (auto const& term : op.tensors.terms) {
        if (term.coefficient == Complex{}) continue;
        for (std::size_t q = 0; q < d.targets.size(); ++q) {
            double total = 0.;
            for (auto const& factor : term.factors) {
                int const site = factor.site, dim = d.dims[static_cast<std::size_t>(site)];
                auto const& matrix = factor.matrix;
                auto const& values = d.spaces[static_cast<std::size_t>(site)].quantum_numbers[q];
                std::optional<double> delta;
                for (int r = 0; r < dim; ++r)
                    for (int c = 0; c < dim; ++c)
                        if (matrix[static_cast<std::size_t>(r * dim + c)] != Complex{}) {
                            double const next = values[static_cast<std::size_t>(r)] - values[static_cast<std::size_t>(c)];
                            if (delta && next != *delta) return false;
                            delta = next;
                        }
                total += delta.value_or(0.);
            }
            if (total != 0.) return false;
        }
    }
    return true;
}
double difference(Action const& a, Action const& b) {
    double result = 0.;
    for (auto const& [state, value] : a) {
        auto const it = b.find(state);
        result = std::max(result, std::abs(value - (it == b.end() ? Complex{} : it->second)));
    }
    for (auto const& [state, value] : b)
        if (a.find(state) == a.end()) result = std::max(result, std::abs(value));
    return result;
}
void validate(Operator const& op, ProjectedData const& d) {
    if (certified(op, d)) return;
    double scale = 0., hermitian = 0., translation_error = 0.;
    bool charge_translation_mismatch = false;
    detail::visit_product_basis(
        d.dims, d.constraints, std::numeric_limits<std::size_t>::max(), [&](State const& state) {
            auto const actual = action(op, d, state, false, true), dagger = action(op, d, state, true);
            auto const source_labels = charge_labels(d, charges(d, state));
            for (auto const& entry : actual)
                scale = std::max({scale, std::abs(entry.second.real()), std::abs(entry.second.imag())});
            hermitian = std::max(hermitian, difference(actual, dagger));
            for (auto const& t : d.translations) {
                State shifted;
                double const input_sign = translate(d, state, t.sites, shifted);
                charge_translation_mismatch =
                    charge_translation_mismatch || charge_labels(d, charges(d, shifted)) != source_labels;
                auto const shifted_action = action(op, d, shifted);
                Action expected;
                for (auto const& entry : actual) {
                    double const output_sign = translate(d, entry.first, t.sites, shifted);
                    expected[shifted] += input_sign * output_sign * entry.second;
                }
                translation_error = std::max(translation_error, difference(shifted_action, expected));
            }
        });
    double const tolerance = 1e-11 * scale;
    if (hermitian > tolerance) invalid("Hamiltonian must be Hermitian");
    if (charge_translation_mismatch) invalid("translation changes a conserved charge cluster");
    if (translation_error > tolerance) invalid("Hamiltonian is not invariant under the lattice translations");
}
struct Transition {
    std::size_t row, col, group;
    Complex value;
};
std::vector<Transition> transitions(Operator const& op, ProjectedData const& d) {
    if (d.sizes.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        invalid("block exceeds sparse index range");
    std::vector<std::vector<std::size_t>> labels;
    labels.reserve(d.sizes.size());
    for (std::size_t i = 0; i < d.sizes.size(); ++i) labels.push_back(charge_labels(d, charges(d, row(d, i))));
    std::vector<Transition> entries;
    for (std::size_t col = 0; col < d.sizes.size(); ++col) {
        auto const source = row(d, col);
        act(op, source, [&](State const& target, Complex value) {
            auto const c = canonical(d, target);
            auto const found = d.indices.find(c.state);
            // Both callers validate the selected parent first. Cross-charge
            // amplitudes are therefore exactly zero, even when termwise double
            // accumulation would leave a residue in a parent mixing charges.
            // Validated translations preserve these representative labels.
            if (found != d.indices.end() && labels[found->second] == labels[col])
                entries.push_back({found->second, col, c.group,
                                   std::sqrt(static_cast<double>(d.sizes[col]) / static_cast<double>(d.sizes[found->second]))
                                       * c.sign * value});
        });
    }
    return entries;
}
Matrix assemble_hamiltonian(std::vector<Transition> const& transitions, ProjectedData const& d,
                            std::vector<int> const* indices = nullptr) {
    CooFromLocal entries;
    for (auto const& t : transitions) {
        auto const r = indices ? (*indices)[t.row] : static_cast<int>(t.row);
        auto const c = indices ? (*indices)[t.col] : static_cast<int>(t.col);
        if (r >= 0 && c >= 0) {
            entries.rows.push_back(r);
            entries.cols.push_back(c);
            entries.data.push_back(t.value * d.characters[t.group]);
        }
    }
    try {
        return canonical_csr<Complex>(entries, d.sizes.size(), d.sizes.size());
    } catch (std::invalid_argument const&) {
        invalid("projected matrix overflowed");
    }
}
template<class Emit>
void operator_entries(Operator const& op, ProjectedData const& source, ProjectedData const& target, Emit&& emit,
                      std::vector<Complex> const* state = nullptr) {
    if (source.dims != target.dims || source.spaces.size() != target.spaces.size())
        invalid("operator bases have different local spaces");
    for (std::size_t site = 0; site < source.spaces.size(); ++site)
        if (source.spaces[site].parity != target.spaces[site].parity
            || source.spaces[site].quantum_numbers != target.spaces[site].quantum_numbers)
            invalid("operator bases have different local spaces");
    for (std::size_t col = 0; col < source.sizes.size(); ++col) {
        if (state && (*state)[col] == Complex{}) continue;
        for (auto const& item : orbit(source, col))
            act(op, item.first, [&](State const& output, Complex value) {
                project_components(target, output, [&](std::size_t row, Complex coefficient) {
                    emit(row, col, std::conj(coefficient) * value * item.second);
                });
            });
    }
}
void check_state(ProjectedData const& d, std::vector<Complex> const& state) {
    if (state.size() != d.sizes.size()) invalid("state length does not match basis dimension");
    for (auto z : state)
        if (!finite(z)) invalid("state amplitudes must be finite");
}
} // namespace

ProjectedBasis::ProjectedBasis(std::shared_ptr<const detail::ProjectedData> data) : data_(std::move(data)) {}
std::size_t ProjectedBasis::dimension() const { return data_->sizes.size(); }
std::size_t ProjectedBasis::num_sites() const { return data_->dims.size(); }
std::size_t ProjectedBasis::parent_dimension() const { return data_->parent_dimension; }
std::size_t ProjectedBasis::momentum_index() const { return data_->momentum_index; }
std::vector<int> const& ProjectedBasis::site_dims() const { return data_->dims; }
std::vector<SpectrumLocalSpace> const& ProjectedBasis::local_spaces() const { return data_->spaces; }
std::vector<SpectrumTranslation> const& ProjectedBasis::translations() const { return data_->translations; }
std::vector<std::optional<double>> const& ProjectedBasis::charge_targets() const { return data_->targets; }
std::vector<double> const& ProjectedBasis::momentum() const { return data_->momentum; }
std::vector<std::int32_t> const& ProjectedBasis::representatives() const { return data_->representatives; }
std::vector<std::size_t> const& ProjectedBasis::orbit_sizes() const { return data_->sizes; }

ProjectedBasis make_projected_basis(std::vector<SpectrumLocalSpace> const& spaces,
                                    std::vector<SpectrumTranslation> const& translations,
                                    std::vector<std::optional<double>> const& targets,
                                    std::size_t momentum_index) {
    auto data = metadata(spaces, translations, targets, momentum_index);
    fill_basis(*data);
    return ProjectedAccess::make(std::move(data));
}

CooFromLocal projected_operator(models::CompiledOperator const& input, ProjectedBasis const& source,
                                ProjectedBasis const& target) {
    auto const& s = ProjectedAccess::get(source);
    auto const& t = ProjectedAccess::get(target);
    auto const op = prepare_operator(input, s);
    CooFromLocal result;
    operator_entries(op, s, t, [&](std::size_t r, std::size_t c, Complex z) {
        result.rows.push_back(static_cast<std::int64_t>(r));
        result.cols.push_back(static_cast<std::int64_t>(c));
        result.data.push_back(z);
        result.is_complex = result.is_complex || z.imag() != 0.;
    });
    return result;
}
std::vector<Complex> apply_projected_operator(models::CompiledOperator const& input, ProjectedBasis const& source,
                                              ProjectedBasis const& target, std::vector<Complex> const& state) {
    auto const& s = ProjectedAccess::get(source);
    auto const& t = ProjectedAccess::get(target);
    check_state(s, state);
    auto const op = prepare_operator(input, s);
    std::vector<Complex> result(t.sizes.size());
    operator_entries(op, s, t, [&](std::size_t r, std::size_t c, Complex z) { result[r] += z * state[c]; }, &state);
    for (auto z : result)
        if (!finite(z)) invalid("projected state action overflowed");
    return result;
}
Complex expect_projected_operator(models::CompiledOperator const& input, ProjectedBasis const& basis,
                                  std::vector<Complex> const& state) {
    auto const output = apply_projected_operator(input, basis, basis, state);
    Complex result = 0.;
    for (std::size_t i = 0; i < state.size(); ++i) result += std::conj(state[i]) * output[i];
    if (!finite(result)) invalid("expectation value overflowed");
    return result;
}
ProjectedExpansion expand_projected_state(ProjectedBasis const& basis, std::vector<Complex> const& state) {
    auto const& d = ProjectedAccess::get(basis);
    check_state(d, state);
    ProjectedExpansion result;
    result.basis.dimension = detail::visit_product_basis(
        d.dims, d.constraints, std::numeric_limits<std::size_t>::max(), [&](State const& product) {
            result.basis.states.insert(result.basis.states.end(), product.begin(), product.end());
            Complex amplitude{};
            project_components(d, product, [&](std::size_t column, Complex coefficient) {
                amplitude += coefficient * state[column];
            });
            result.amplitudes.push_back(amplitude);
        });
    return result;
}

struct ProjectedProblem::Data {
    ProjectedBasis basis;
    std::vector<std::string> charges;
    Matrix matrix;
};
ProjectedProblem::ProjectedProblem(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
ProjectedProblem::Data const& ProjectedProblem::data() const {
    if (!data_) throw std::logic_error("ProjectedProblem: moved-from problem");
    return *data_;
}
ProjectedBasis const& ProjectedProblem::basis() const { return data().basis; }
std::size_t ProjectedProblem::dimension() const { return data().basis.dimension(); }
std::vector<std::string> const& ProjectedProblem::charges() const { return data().charges; }
SparseMatrix<Complex> const& ProjectedProblem::hamiltonian() const { return data().matrix; }
ProjectedProblem ProjectedProblem::prepare(models::CompiledOperator const& input, ProjectedBasis const& basis,
                                           std::vector<std::string> charge_names) {
    auto const& d = ProjectedAccess::get(basis);
    auto const op = prepare_operator(input, d);
    validate(op, d);
    auto matrix = assemble_hamiltonian(transitions(op, d), d);
    return ProjectedProblem(std::make_shared<const Data>(Data{basis, std::move(charge_names), std::move(matrix)}));
}
ProjectedProblem prepare_projected(models::CompiledOperator const& input, ProjectedBasis const& basis) {
    return ProjectedProblem::prepare(input, basis, {});
}
ProjectedEigenResult ProjectedProblem::eigenstates(EigenstateOptions const& options) const {
    auto const& source = data();
    detail::validate_request(options);
    auto const plan = detail::eigenstate_plan(options, dimension());
    auto solved = detail::solve_block(source.matrix, options.which, plan, true);
    solved.spectrum.momentum = source.basis.momentum();
    return {*this, std::move(solved.spectrum), std::move(solved.vectors), std::move(solved.residuals),
            detail::lanczos_requested(options.method, dimension())};
}

std::vector<SpectrumBlock> symmetry_spectrum_from_local(
    models::CompiledOperator const& input, std::vector<SpectrumLocalSpace> const& spaces,
    std::vector<SpectrumTranslation> const& translations, std::vector<std::optional<double>> const& targets,
    int count, std::vector<std::size_t> const& requested, std::string const& which,
    std::optional<double> tolerance, std::optional<int> max_iterations) {
    validate_eigensolver_controls(which, tolerance, max_iterations);
    if (count <= 0) invalid("eigenvalue count must be positive");
    auto prototype = metadata(spaces, translations, targets, 0);
    auto const op = prepare_operator(input, *prototype);
    std::vector<std::size_t> momenta(requested);
    if (momenta.empty()) {
        momenta.resize(prototype->permutations.size());
        std::iota(momenta.begin(), momenta.end(), std::size_t{0});
    } else {
        std::sort(momenta.begin(), momenta.end());
        if (momenta.back() >= prototype->permutations.size()
            || std::adjacent_find(momenta.begin(), momenta.end()) != momenta.end())
            invalid("requested momenta must be unique and inside translation grid");
    }
    fill_basis(*prototype, true);
    validate(op, *prototype);
    auto const orbit_transitions = transitions(op, *prototype);
    std::vector<SpectrumBlock> result;
    std::size_t projected_dimension = 0;
    for (auto k : momenta) {
        auto d = metadata(spaces, translations, targets, k);
        d->parent_dimension = prototype->parent_dimension;
        d->charge_anchors = prototype->charge_anchors;
        std::vector<int> orbit_indices(prototype->sizes.size(), -1);
        State scratch;
        for (std::size_t old = 0; old < prototype->sizes.size(); ++old) {
            auto const rep = row(*prototype, old);
            Complex sum = 0.;
            for (std::size_t g = 0; g < d->permutations.size(); ++g) {
                double const sign = translate(*d, rep, d->permutations[g], scratch);
                if (scratch == rep) sum += sign * d->characters[g];
            }
            if (std::abs(sum) <= 64. * std::numeric_limits<double>::epsilon() * static_cast<double>(d->permutations.size()))
                continue;
            auto const index = d->sizes.size();
            orbit_indices[old] = static_cast<int>(index);
            d->sizes.push_back(prototype->sizes[old]);
            d->indices.emplace(rep, index);
            d->representatives.insert(d->representatives.end(), rep.begin(), rep.end());
        }
        projected_dimension += d->sizes.size();
        auto const matrix = assemble_hamiltonian(orbit_transitions, *d, &orbit_indices);
        std::map<std::vector<std::size_t>, std::vector<std::size_t>> groups;
        for (std::size_t r = 0; r < d->sizes.size(); ++r) groups[charge_labels(*d, charges(*d, row(*d, r)))].push_back(r);
        for (auto const& group : groups) {
            auto const& label = group.first;
            auto const& members = group.second;
            std::vector<int> indices(d->sizes.size(), -1);
            for (std::size_t i = 0; i < members.size(); ++i) indices[members[i]] = static_cast<int>(i);
            CooFromLocal entries;
            for (auto r : members)
                for (auto p = matrix.offsets[r]; p < matrix.offsets[r + 1]; ++p) {
                    auto const column = static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(p)]);
                    if (indices[column] < 0) continue;
                    entries.rows.push_back(indices[r]);
                    entries.cols.push_back(indices[column]);
                    entries.data.push_back(matrix.values[static_cast<std::size_t>(p)]);
                }
            auto const block = canonical_csr<Complex>(entries, members.size(), members.size());
            auto solved = detail::solve_block(block, which,
                detail::spectrum_plan(members.size(), count, which, tolerance, max_iterations), false);
            solved.spectrum.momentum = d->momentum;
            for (std::size_t q = 0; q < label.size(); ++q)
                solved.spectrum.quantum_numbers.push_back(d->charge_anchors[q][label[q]]);
            result.push_back(std::move(solved.spectrum));
        }
    }
    if (momenta.size() == prototype->permutations.size() && projected_dimension != prototype->parent_dimension)
        throw std::runtime_error("projected ED: momentum bases do not span the selected product sector");
    return result;
}
} // namespace alps::ed
