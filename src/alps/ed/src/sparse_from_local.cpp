// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/assembly/sparse_from_local.cpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names.
//
// Implementation of alps::ed::build_sparse_from_local (see the header for
// the locked conventions).

#include <alps/ed/sparse_from_local.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include "local_action.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

namespace alps {
namespace ed {

namespace {
using detail::PreparedTerm;
using detail::prepare_terms;
using detail::apply_term;

double real_or_throw(std::complex<double> const& z, char const* what) {
    if (std::abs(z.imag()) > 1e-12)
        throw std::invalid_argument(
            std::string("build_sparse_from_local: complex ") + what +
            " — this is the real CSR path");
    return z.real();
}

// The real CSR path uses site 0 least significant; the complex visitor uses
// the opposite strides. Neither representation determines matrix order:
// code lookup retains exactly the caller's row ordering.
struct EncodedBasis {
    using State = std::uint64_t;
    std::vector<int> const& dims;
    std::vector<State> const& strides;
    std::vector<State> const& codes;
    std::unordered_map<State, std::int64_t> const* indices = nullptr;
    std::int64_t const* local_states = nullptr;   // row-major rows, or null

    State state(std::size_t row) const { return codes[row]; }
    int digit(std::size_t row, int site) const {
        auto const s = static_cast<std::size_t>(site);
        if (local_states) return static_cast<int>(local_states[row * dims.size() + s]);
        return static_cast<int>((codes[row] / strides[s]) % static_cast<State>(dims[s]));
    }
    bool unchanged(State target, std::size_t input) const { return target == codes[input]; }
    void shift(State& target, int site, int before, int after) const {
        // Unsigned arithmetic avoids signed overflow near the code limit.
        target -= static_cast<State>(before) * strides[static_cast<std::size_t>(site)];
        target += static_cast<State>(after) * strides[static_cast<std::size_t>(site)];
    }
    std::int64_t index(State target) const {
        if (indices) {
            auto const found = indices->find(target);
            return found == indices->end() ? -1 : found->second;
        }
        auto const found = std::lower_bound(codes.begin(), codes.end(), target);
        return found == codes.end() || *found != target ? -1 : found - codes.begin();
    }
};

// When the full product does not fit uint64, compare explicit rows. This
// still only stores and searches the selected basis, never the full space.
struct ExplicitBasis {
    using State = std::vector<std::int64_t>;
    std::int64_t const* values;
    std::size_t width;
    std::vector<std::size_t> order;

    std::int64_t const* row(std::size_t i) const { return values + i * width; }
    State state(std::size_t i) const {
        auto const r = row(i);
        return State(r, r + width);
    }
    int digit(std::size_t i, int site) const { return static_cast<int>(row(i)[site]); }
    void shift(State& target, int site, int, int after) const { target[static_cast<std::size_t>(site)] = after; }
    bool unchanged(State const& target, std::size_t input) const {
        auto const r = row(input);
        return std::equal(r, r + width, target.begin());
    }
    std::int64_t index(State const& target) const {
        auto const found = std::lower_bound(order.begin(), order.end(), target,
            [&](std::size_t i, State const& other) {
                auto const r = row(i);
                return std::lexicographical_compare(r, r + width, other.begin(), other.end());
            });
        if (found == order.end()) return -1;
        auto const r = row(*found);
        return std::equal(r, r + width, target.begin()) ? static_cast<std::int64_t>(*found) : -1;
    }
};

template <typename Basis>
bool visit_terms(std::vector<PreparedTerm<std::complex<double>>> const& terms,
                 Basis const& basis, std::size_t count,
                 std::vector<std::vector<double> const*> const& parity,
                 MatrixElementVisitor const& visitor) {
    bool is_complex = false;
    for (auto const& term : terms) {
        if (term.coeff == std::complex<double>{}) continue;
        for (std::size_t col = 0; col < count; ++col)
            apply_term(term, basis, parity, col, [&](std::int64_t row, std::complex<double> value) {
                if (value == std::complex<double>{}) return;
                visitor(row, static_cast<std::int64_t>(col), value);
                is_complex = is_complex || term.coeff.imag() != 0. || value.imag() != 0.;
            });
    }
    return is_complex;
}

struct OperatorShape {
    std::vector<std::uint64_t> strides;
    std::uint64_t product = 1;
    bool encodable = true;
    std::vector<std::vector<double> const*> parity;
};

OperatorShape prepare_shape(models::CompiledOperator const& input,
                            std::vector<std::vector<double>> const& parity_values) {
    auto const n = input.site_dims.size();
    for (int d : input.site_dims)
        if (d < 1) throw std::invalid_argument("sparse assembly: site dim < 1");
    if (!parity_values.empty() && parity_values.size() != n)
        throw std::invalid_argument("sparse assembly: parity site count mismatch");
    OperatorShape shape;
    shape.parity.assign(n, nullptr);
    for (std::size_t site = 0; site < parity_values.size(); ++site) {
        if (parity_values[site].size() != static_cast<std::size_t>(input.site_dims[site]))
            throw std::invalid_argument("sparse assembly: parity length != site dim");
        shape.parity[site] = &parity_values[site];
    }
    shape.strides.resize(n);
    for (std::size_t site = n; site-- > 0;) {
        shape.strides[site] = shape.product;
        auto const d = static_cast<std::uint64_t>(input.site_dims[site]);
        if (shape.product > std::numeric_limits<std::uint64_t>::max() / d) {
            shape.encodable = false;
            break;
        }
        shape.product *= d;
    }
    return shape;
}

struct CodeIndex {
    std::unordered_map<std::uint64_t, std::int64_t> values;
    bool sorted;

    explicit CodeIndex(std::vector<std::uint64_t> const& codes)
        : sorted(std::is_sorted(codes.begin(), codes.end())) {
        if (sorted) {
            if (std::adjacent_find(codes.begin(), codes.end()) != codes.end())
                throw std::invalid_argument("sparse assembly: basis states must be distinct");
        } else {
            values.reserve(codes.size());
            for (std::size_t row = 0; row < codes.size(); ++row)
                if (!values.emplace(codes[row], static_cast<std::int64_t>(row)).second)
                    throw std::invalid_argument("sparse assembly: basis states must be distinct");
        }
    }

    auto const* lookup() const { return sorted ? nullptr : &values; }
};

} // namespace

bool visit_local_transitions(
    models::CompiledOperator const& input,
    std::vector<std::int64_t> const& sources, std::size_t num_sources,
    LocalTransitionVisitor const& visitor,
    std::vector<std::vector<double>> const& parity_values) {
    if (!visitor) throw std::invalid_argument("sparse assembly: empty local transition visitor");
    auto const n = input.site_dims.size();
    if (num_sources > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()) ||
        (n && num_sources > std::numeric_limits<std::size_t>::max() / n) ||
        sources.size() != num_sources * n)
        throw std::invalid_argument("sparse assembly: source shape does not match site count");
    auto const shape = prepare_shape(input, parity_values);
    for (std::size_t row = 0; row < num_sources; ++row)
        for (std::size_t site = 0; site < n; ++site) {
            auto const local = sources[row * n + site];
            if (local < 0 || local >= input.site_dims[site])
                throw std::invalid_argument("sparse assembly: source local state out of range");
        }
    auto const terms = prepare_terms<std::complex<double>>(
        input, [](auto value, char const*) { return value; }, true);
    ExplicitBasis const basis{sources.data(), n, {}};
    bool is_complex = false;
    for (auto const& term : terms) {
        if (term.coeff == std::complex<double>{}) continue;
        for (std::size_t source = 0; source < num_sources; ++source)
            detail::apply_term_state(term, basis, shape.parity, source,
                [&](auto const& target, std::complex<double> amplitude) {
                    if (amplitude == std::complex<double>{}) return;
                    visitor(static_cast<std::int64_t>(source), target, amplitude);
                    is_complex = is_complex || term.coeff.imag() != 0. || amplitude.imag() != 0.;
                });
    }
    return is_complex;
}

bool visit_local_operator(
    models::CompiledOperator const& input,
    std::vector<std::int64_t> const& states, std::size_t num_states,
    MatrixElementVisitor const& visitor,
    std::vector<std::vector<double>> const& parity_values) {
    if (!visitor) throw std::invalid_argument("sparse assembly: empty matrix element visitor");
    auto const n = input.site_dims.size();
    if (num_states > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()) ||
        (n && num_states > std::numeric_limits<std::size_t>::max() / n) ||
        states.size() != num_states * n)
        throw std::invalid_argument("sparse assembly: basis shape does not match site count");
    auto const shape = prepare_shape(input, parity_values);
    for (std::size_t row = 0; row < num_states; ++row)
        for (std::size_t site = 0; site < n; ++site) {
            auto const local = states[row * n + site];
            if (local < 0 || local >= input.site_dims[site])
                throw std::invalid_argument("sparse assembly: basis local state out of range");
        }
    auto const terms = prepare_terms<std::complex<double>>(
        input, [](auto value, char const*) { return value; }, true);

    // A fast code index when possible; otherwise a lexicographic index over
    // explicit rows. Both indices return the caller's original row number.
    if (shape.encodable) {
        std::vector<std::uint64_t> codes(num_states, 0);
        for (std::size_t row = 0; row < num_states; ++row)
            for (std::size_t site = 0; site < n; ++site)
                codes[row] += static_cast<std::uint64_t>(states[row * n + site]) * shape.strides[site];
        CodeIndex const indices(codes);
        return visit_terms(terms,
                           EncodedBasis{input.site_dims, shape.strides, codes, indices.lookup(), states.data()},
                           num_states, shape.parity, visitor);
    }
    ExplicitBasis basis{states.data(), n, std::vector<std::size_t>(num_states)};
    std::iota(basis.order.begin(), basis.order.end(), 0);
    std::sort(basis.order.begin(), basis.order.end(), [&](auto i, auto j) {
        auto const a = basis.row(i), b = basis.row(j);
        return std::lexicographical_compare(a, a + n, b, b + n);
    });
    for (std::size_t i = 1; i < num_states; ++i) {
        auto const a = basis.row(basis.order[i - 1]), b = basis.row(basis.order[i]);
        if (std::equal(a, a + n, b))
            throw std::invalid_argument("sparse assembly: basis states must be distinct");
    }
    return visit_terms(terms, basis, num_states, shape.parity, visitor);
}

bool visit_local_operator(
    models::CompiledOperator const& input, std::vector<std::uint64_t> const& codes,
    MatrixElementVisitor const& visitor,
    std::vector<std::vector<double>> const& parity_values) {
    if (!visitor) throw std::invalid_argument("sparse assembly: empty matrix element visitor");
    if (codes.size() > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()))
        throw std::invalid_argument("sparse assembly: basis has too many states");
    auto const shape = prepare_shape(input, parity_values);
    if (!shape.encodable)
        throw std::invalid_argument("sparse assembly: full Hilbert space does not fit uint64_t codes");
    for (auto code : codes)
        if (code >= shape.product)
            throw std::invalid_argument("sparse assembly: basis code out of range");
    auto const terms = prepare_terms<std::complex<double>>(
        input, [](auto value, char const*) { return value; }, true);
    CodeIndex const indices(codes);
    return visit_terms(terms,
                       EncodedBasis{input.site_dims, shape.strides, codes, indices.lookup(), nullptr},
                       codes.size(), shape.parity, visitor);
}

CooFromLocal build_coo_from_local(
    models::CompiledOperator const& input,
    std::vector<std::int64_t> const& states, std::size_t num_states,
    std::vector<std::vector<double>> const& parity_values) {
    CooFromLocal out;
    out.is_complex = visit_local_operator(input, states, num_states,
        [&](std::int64_t row, std::int64_t col, std::complex<double> value) {
            out.rows.push_back(row);
            out.cols.push_back(col);
            out.data.push_back(value);
        }, parity_values);
    return out;
}

namespace {

// Row r of A is the conjugate of column r of A^dagger, so applying the adjoint
// to one basis state yields one output row of A. The adjoint conjugates each
// coefficient and conjugate-transposes each factor. A string sign F_s(in_s) on
// a factor site depends on A's input, which is the adjoint's output: it moves
// into the factor's rows. Strings on other sites see an unchanged digit and
// stay. Conjugation and +-1 factors are exact, so every contribution is the
// exact conjugate of the one the column-wise kernel emits.
models::CompiledOperator row_operator(models::CompiledOperator const& input,
                                      std::vector<std::vector<double>> const& parity) {
    auto const n = input.site_dims.size();
    models::CompiledOperator out;
    out.site_dims = input.site_dims;
    out.terms.reserve(input.terms.size());
    for (auto const& term : input.terms) {
        models::LocalTerm adjoint;
        adjoint.coefficient = std::conj(term.coefficient);
        std::vector<int> folded;
        for (auto const& factor : term.factors) {
            if (factor.site < 0 || static_cast<std::size_t>(factor.site) >= n)
                throw std::invalid_argument("sparse assembly: factor site out of range");
            auto const d = static_cast<std::size_t>(input.site_dims[static_cast<std::size_t>(factor.site)]);
            if (d > std::numeric_limits<std::size_t>::max() / d || factor.matrix.size() != d * d)
                throw std::invalid_argument("sparse assembly: operator size != site dim^2");
            auto const strings = parity.empty() ? 0 : std::count(term.jordan_wigner.begin(),
                                                                 term.jordan_wigner.end(), factor.site);
            if (strings) folded.push_back(factor.site);
            models::LocalFactor dagger{factor.site, std::vector<std::complex<double>>(d * d)};
            for (std::size_t a = 0; a < d; ++a)
                for (std::size_t b = 0; b < d; ++b) {
                    auto value = std::conj(factor.matrix[b * d + a]);
                    if (strings % 2) value *= parity[static_cast<std::size_t>(factor.site)][a];
                    dagger.matrix[a * d + b] = value;
                }
            adjoint.factors.push_back(std::move(dagger));
        }
        for (int site : term.jordan_wigner)
            if (std::find(folded.begin(), folded.end(), site) == folded.end()) adjoint.jordan_wigner.push_back(site);
        out.terms.push_back(std::move(adjoint));
    }
    return out;
}

bool has_complex_entries(models::CompiledOperator const& input) {
    for (auto const& term : input.terms) {
        if (term.coefficient.imag() != 0.) return true;
        for (auto const& factor : term.factors)
            for (auto value : factor.matrix)
                if (value.imag() != 0.) return true;
    }
    return false;
}

// One row's contributions sorted stably by column (term order within a
// column, as canonical_csr's stable sort keeps COO order), then reduced.
template <typename Scalar>
void append_row(std::vector<std::pair<std::int64_t, std::complex<double>>>& entries, SparseMatrix<Scalar>& result) {
    std::stable_sort(entries.begin(), entries.end(), [](auto const& a, auto const& b) { return a.first < b.first; });
    auto first = entries.begin();
    while (first != entries.end()) {
        auto const column = first->first;
        Scalar value = detail::matrix_cast<Scalar>(first++->second);
        while (first != entries.end() && first->first == column)
            value = detail::matrix_add(value, detail::matrix_cast<Scalar>(first++->second));
        if (!detail::finite_scalar(value)) throw std::invalid_argument("sparse matrix accumulation overflowed");
        if (value != Scalar{}) { result.indices.push_back(column); result.values.push_back(value); }
    }
    result.offsets.push_back(static_cast<std::int64_t>(result.values.size()));
}

template <typename Scalar, typename Basis>
bool assemble_rows(std::vector<PreparedTerm<std::complex<double>>> const& terms, Basis const& basis,
                   std::size_t count, std::vector<std::vector<double> const*> const& parity,
                   SparseMatrix<Scalar>& result) {
    bool is_complex = false;
    result.rows = result.columns = count;
    result.offsets.reserve(count + 1);
    result.offsets.push_back(0);
    std::vector<std::pair<std::int64_t, std::complex<double>>> entries;
    for (std::size_t row = 0; row < count; ++row) {
        entries.clear();
        for (auto const& term : terms) {
            if (term.coeff == std::complex<double>{}) continue;
            apply_term(term, basis, parity, row, [&](std::int64_t column, std::complex<double> value) {
                if (value == std::complex<double>{}) return;
                if (!detail::finite_scalar(value)) throw std::invalid_argument("sparse matrix contribution must be finite");
                entries.emplace_back(column, std::conj(value));
                is_complex = is_complex || term.coeff.imag() != 0. || value.imag() != 0.;
            });
        }
        append_row(entries, result);
    }
    return is_complex;
}

template <typename Scalar>
bool assemble_operator(models::CompiledOperator const& adjoint, std::vector<std::int64_t> const& states,
                       std::size_t count, OperatorShape const& shape, SparseMatrix<Scalar>& result) {
    auto const n = adjoint.site_dims.size();
    auto const terms = prepare_terms<std::complex<double>>(
        adjoint, [](auto value, char const*) { return value; }, true);
    if (shape.encodable) {
        std::vector<std::uint64_t> codes(count, 0);
        for (std::size_t row = 0; row < count; ++row)
            for (std::size_t site = 0; site < n; ++site)
                codes[row] += static_cast<std::uint64_t>(states[row * n + site]) * shape.strides[site];
        CodeIndex const indices(codes);
        return assemble_rows(terms, EncodedBasis{adjoint.site_dims, shape.strides, codes, indices.lookup(), states.data()},
                             count, shape.parity, result);
    }
    ExplicitBasis basis{states.data(), n, std::vector<std::size_t>(count)};
    std::iota(basis.order.begin(), basis.order.end(), 0);
    std::sort(basis.order.begin(), basis.order.end(), [&](auto i, auto j) {
        auto const a = basis.row(i), b = basis.row(j);
        return std::lexicographical_compare(a, a + n, b, b + n);
    });
    for (std::size_t i = 1; i < count; ++i) {
        auto const a = basis.row(basis.order[i - 1]), b = basis.row(basis.order[i]);
        if (std::equal(a, a + n, b)) throw std::invalid_argument("sparse assembly: basis states must be distinct");
    }
    return assemble_rows(terms, basis, count, shape.parity, result);
}

} // namespace

CanonicalCsr build_canonical_csr(models::CompiledOperator const& input,
    std::vector<std::int64_t> const& states, std::size_t count,
    std::vector<std::vector<double>> const& parity_values, bool complex_values) {
    auto const n = input.site_dims.size();
    if (count >= static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max()) ||
        (n && count > std::numeric_limits<std::size_t>::max() / n) || states.size() != count * n)
        throw std::invalid_argument("sparse assembly: basis shape does not match site count");
    auto const shape = prepare_shape(input, parity_values);
    for (std::size_t row = 0; row < count; ++row)
        for (std::size_t site = 0; site < n; ++site) {
            auto const local = states[row * n + site];
            if (local < 0 || local >= input.site_dims[site])
                throw std::invalid_argument("sparse assembly: basis local state out of range");
        }
    auto const adjoint = row_operator(input, parity_values);
    // A real operator emits only real contributions: hold real values from the start.
    if (!complex_values && !has_complex_entries(input)) {
        SparseMatrix<double> result;
        assemble_operator(adjoint, states, count, shape, result);
        return result;
    }
    SparseMatrix<std::complex<double>> result;
    bool const is_complex = assemble_operator(adjoint, states, count, shape, result);
    if (complex_values || is_complex) return result;
    // Complex entries that never contributed (outside the basis, or zero
    // products): the COO path calls this real, so return the real parts.
    SparseMatrix<double> real{result.rows, result.columns, std::move(result.offsets), std::move(result.indices), {}};
    real.values.reserve(result.values.size());
    for (auto value : result.values) real.values.push_back(value.real());
    return real;
}

SparseFromLocal build_sparse_from_local(
    models::CompiledModel const& model,
    std::map<int, std::vector<double>> const& parity_by_type,
    std::vector<QnConstraint> const& constraints, std::uint64_t max_full_dim) {
    auto const &input = model.op;
    std::size_t const n = input.site_dims.size();
    if (n == 0) throw std::invalid_argument("build_sparse_from_local: no sites");
    if (model.site_types.size() != n)
        throw std::invalid_argument(
            "build_sparse_from_local: site_types/site_dims length mismatch");

    // Mixed-radix layout, site 0 least significant.
    std::vector<std::uint64_t> strides(n, 1);
    std::uint64_t full_dim = 1;
    for (std::size_t s = 0; s < n; ++s) {
        int const d = input.site_dims[s];
        if (d < 1)
            throw std::invalid_argument("build_sparse_from_local: site dim < 1");
        strides[s] = full_dim;
        if (full_dim > max_full_dim / static_cast<std::uint64_t>(d))
            throw std::invalid_argument(
                "build_sparse_from_local: full Hilbert space exceeds the cap");
        full_dim *= static_cast<std::uint64_t>(d);
    }
    auto const digit = [&](std::uint64_t state, std::size_t site) -> int {
        return static_cast<int>((state / strides[site]) %
                                static_cast<std::uint64_t>(input.site_dims[site]));
    };

    // Per-site per-state QN values for each constraint (a site type without
    // the QN contributes the constant 0 — the kit convention), and the
    // per-site diagonal parity for the Jordan-Wigner strings.
    std::vector<std::vector<std::vector<double>>> qn_values(constraints.size());
    for (std::size_t c = 0; c < constraints.size(); ++c) {
        qn_values[c].resize(n);
        for (std::size_t s = 0; s < n; ++s) {
            auto const ti = model.diag_observables.find(model.site_types[s]);
            if (ti != model.diag_observables.end()) {
                auto const qi = ti->second.find(constraints[c].qn);
                if (qi != ti->second.end()) {
                    if (static_cast<int>(qi->second.size()) != input.site_dims[s])
                        throw std::invalid_argument(
                            "build_sparse_from_local: QN '" + constraints[c].qn +
                            "' value count != site dim");
                    qn_values[c][s] = qi->second;
                    continue;
                }
            }
            qn_values[c][s].assign(static_cast<std::size_t>(input.site_dims[s]),
                                   0.0);
        }
    }
    std::vector<std::vector<double> const*> parity(n, nullptr);
    for (std::size_t s = 0; s < n; ++s) {
        auto const pi = parity_by_type.find(model.site_types[s]);
        if (pi != parity_by_type.end()) {
            if (static_cast<int>(pi->second.size()) != input.site_dims[s])
                throw std::invalid_argument(
                    "build_sparse_from_local: parity length != site dim");
            parity[s] = &pi->second;
        }
    }

    // Basis enumeration with the sector filter.
    SparseFromLocal out;
    for (std::uint64_t state = 0; state < full_dim; ++state) {
        bool keep = true;
        for (std::size_t c = 0; c < constraints.size() && keep; ++c) {
            double sum = 0.0;
            for (std::size_t s = 0; s < n; ++s)
                sum += qn_values[c][s][static_cast<std::size_t>(digit(state, s))];
            keep = std::abs(sum - constraints[c].target) < 1e-9;
        }
        if (keep) out.states.push_back(state);
    }
    out.dim = static_cast<std::int64_t>(out.states.size());

    auto const terms = prepare_terms<double>(input, real_or_throw, false);
    EncodedBasis const basis{input.site_dims, strides, out.states, nullptr, nullptr};

    // Row-by-row scatter.
    out.indptr.reserve(static_cast<std::size_t>(out.dim) + 1);
    out.indptr.push_back(0);
    std::map<std::int64_t, double> row;
    for (std::int64_t r = 0; r < out.dim; ++r) {
        row.clear();
        for (auto const& term : terms)
            apply_term(term, basis, parity, static_cast<std::size_t>(r),
                       [&](std::int64_t col, double value) { row[col] += value; });
        for (auto const& [col, v] : row) {
            out.indices.push_back(col);
            out.data.push_back(v);
        }
        out.indptr.push_back(static_cast<std::int64_t>(out.indices.size()));
    }
    return out;
}

} // namespace ed
} // namespace alps
