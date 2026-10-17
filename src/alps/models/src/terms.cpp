// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/operators/terms.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
//
// Shared native operator algebra. Ordinary-model conventions are pinned by
// models_kit_parity.cpp; term_accumulation.cpp locks numerical preparation.

#include <alps/models/terms.hpp>
#include "exact_arithmetic.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <unordered_map>

namespace alps {
namespace models {

namespace {
using detail::ExactComplex;
using ExactMatrix = std::vector<ExactComplex>;

bool matrix_is_identity(Matrix const& m, int dim, double tol = 1e-12) {
    for (int r = 0; r < dim; ++r)
        for (int c = 0; c < dim; ++c) {
            Complex const expect = (r == c) ? Complex(1.0) : Complex(0.0);
            if (std::abs(m[static_cast<std::size_t>(r) * dim + c] - expect) > tol)
                return false;
        }
    return true;
}

ExactMatrix exact_matrix(Matrix const& matrix) {
    ExactMatrix result;
    result.reserve(matrix.size());
    for (auto value : matrix) result.emplace_back(value);
    return result;
}

ExactMatrix multiply_run(std::vector<Factor> const& factors, std::size_t first,
                         std::size_t last, int dim) {
    auto product = exact_matrix(factors[first].op.matrix);
    for (auto i = first + 1; i < last; ++i) {
        auto const right = exact_matrix(factors[i].op.matrix);
        ExactMatrix next(product.size());
        for (int r = 0; r < dim; ++r) for (int k = 0; k < dim; ++k) {
            auto const& left = product[std::size_t(r) * dim + k];
            if (left.zero()) continue;
            for (int c = 0; c < dim; ++c) {
                auto const& value = right[std::size_t(k) * dim + c];
                if (!value.zero()) next[std::size_t(r) * dim + c] += left * value;
            }
        }
        product = std::move(next);
    }
    return product;
}

bool export_matrix(ExactMatrix const& input, Matrix& output, std::int64_t exponent = 0) {
    output.resize(input.size());
    for (std::size_t i = 0; i < input.size(); ++i)
        if (!detail::export_complex(detail::scale_power_two(input[i], -exponent), output[i]))
            return false;
    return true;
}

// A term-local gauge only: do not materialize a global tensor product. Single
// factors participate when a computed run needs a scale recipient, but the
// ordinary path does not allocate exact matrices for them.
Complex export_scaled_products(std::vector<Factor>& factors,
                                std::vector<ExactMatrix>& products, Complex coefficient) {
    std::vector<detail::ScaleInterval> intervals(factors.size());
    std::vector<std::int64_t> exponents(factors.size());
    std::int64_t sum = 0;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        if (products[i].empty()) products[i] = exact_matrix(factors[i].op.matrix);
        for (auto const& value : products[i]) intervals[i].include(value);
        auto const& bounds = intervals[i];
        if (bounds.lower > bounds.upper)
            throw std::overflow_error("local product components have no finite preserving scale");
        exponents[i] = std::clamp(std::int64_t(0), bounds.lower, bounds.upper);
        sum = detail::checked_exponent_sum(sum, exponents[i]);
    }
    ExactComplex const exact_coefficient(coefficient);
    detail::ScaleInterval coefficient_bounds;
    coefficient_bounds.include(exact_coefficient);
    // Matrices use 2^-e, whereas their common coefficient uses 2^sum(e).
    auto const target = std::clamp(sum, -coefficient_bounds.upper, -coefficient_bounds.lower);
    if (sum == std::numeric_limits<std::int64_t>::lowest())
        throw std::overflow_error("local product scale adjustment exceeds the supported range");
    auto remaining = detail::checked_exponent_sum(target, -sum);
    if (remaining == std::numeric_limits<std::int64_t>::lowest())
        throw std::overflow_error("local product scale adjustment exceeds the supported range");
    for (std::size_t i = 0; i < factors.size() && remaining; ++i) {
        if (remaining > 0) {
            auto const amount = std::min(remaining, intervals[i].upper - exponents[i]);
            exponents[i] += amount;
            remaining -= amount;
        } else {
            auto const amount = std::min(-remaining, exponents[i] - intervals[i].lower);
            exponents[i] -= amount;
            remaining += amount;
        }
    }
    Complex result;
    if (remaining || !detail::export_complex(detail::scale_power_two(exact_coefficient, target), result))
        throw std::overflow_error("local product and coefficient have no finite preserving scale");
    for (std::size_t i = 0; i < factors.size(); ++i)
        if (!export_matrix(products[i], factors[i].op.matrix, exponents[i]))
            throw std::overflow_error("local product cannot be exported without range loss");
    return result;
}

// Like-term key: per factor (site, numeric matrix components, parity).
// Conjugation can change a zero's sign without changing the operator.
std::string term_key(Term const& t) {
    std::string key;
    for (Factor const& f : t.factors) {
        key.append(reinterpret_cast<char const*>(&f.site), sizeof(f.site));
        int const dim = f.op.dim;
        key.append(reinterpret_cast<char const*>(&dim), sizeof(dim));
        for (auto value : f.op.matrix) for (double component : {value.real(), value.imag()}) {
            if (component == 0.) component = 0.;
            key.append(reinterpret_cast<char const*>(&component), sizeof(component));
        }
        char const p = f.op.parity == Parity::odd ? 'o' : 'e';
        key.push_back(p);
    }
    return key;
}

} // namespace

std::map<std::string, double> Term::qn_delta() const {
    std::map<std::string, double> total;
    if (factors.empty()) return total;
    // Intersection of QN names every factor respects.
    bool first = true;
    for (Factor const& f : factors) {
        if (first) {
            for (auto const& [k, v] : f.op.qn_deltas) { (void)v; total[k] = 0.0; }
            first = false;
        } else {
            for (auto it = total.begin(); it != total.end();) {
                if (!f.op.qn_deltas.count(it->first)) it = total.erase(it);
                else ++it;
            }
        }
        if (total.empty()) return total;
    }
    for (Factor const& f : factors)
        for (auto& [k, v] : total) v += f.op.qn_deltas.at(k);
    return total;
}

Term canonicalise_term(Term const& t, bool apply_swap_sign) {
    return canonicalise_term(t, apply_swap_sign, 1e-12);
}

Term canonicalise_term(Term const& t, bool apply_swap_sign, double identity_tolerance) {
    if (!std::isfinite(identity_tolerance) || identity_tolerance < 0)
        throw std::invalid_argument("identity tolerance must be finite and nonnegative");
    if (!std::isfinite(t.coeff.real()) || !std::isfinite(t.coeff.imag()))
        throw std::invalid_argument("canonicalization requires finite coefficients");
    for (auto const& factor : t.factors) {
        if (factor.op.dim <= 0 || factor.op.matrix.size() != std::size_t(factor.op.dim) * factor.op.dim)
            throw std::invalid_argument("canonicalization requires matching positive matrix dimensions");
        for (auto value : factor.op.matrix)
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                throw std::invalid_argument("canonicalization requires finite matrix entries");
    }
    if (t.coeff == Complex{}) return Term{Complex{}, {}, t.sign};
    std::vector<Factor> factors = t.factors;
    int sign = t.sign;

    // Stable bubble sort by site index, adjacent swaps, tracking parity —
    // byte-for-byte the Python loop.
    std::size_t const n = factors.size();
    for (std::size_t i = 1; i < n; ++i) {
        std::size_t j = i;
        while (j > 0) {
            Factor& left = factors[j - 1];
            Factor& right = factors[j];
            if (left.site <= right.site) break;
            if (apply_swap_sign && left.op.parity == Parity::odd &&
                right.op.parity == Parity::odd)
                sign = -sign;
            std::swap(factors[j - 1], factors[j]);
            --j;
        }
    }

    // Keep every multiplication in a complete same-site run exact. Only the
    // final matrix is exported; normal double rounding there remains explicit.
    std::vector<Factor> out;
    std::vector<ExactMatrix> products;
    bool needs_scale = false;
    std::size_t i = 0;
    while (i < factors.size()) {
        int const site_idx = factors[i].site;
        SiteBasisPtr const basis = factors[i].basis;
        int const dim = factors[i].op.dim;
        bool odd = factors[i].op.parity == Parity::odd;
        std::string merged_name = factors[i].op.name;
        std::size_t j = i + 1;
        while (j < factors.size() && factors[j].site == site_idx) {
            if (factors[j].op.dim != dim)
                throw std::invalid_argument("same-site operator dimensions differ");
            if (factors[j].op.parity == Parity::odd) odd = !odd;
            merged_name += "*" + factors[j].op.name;
            ++j;
        }
        Factor merged;
        merged.site = site_idx;
        merged.basis = basis;
        merged.op.name = std::move(merged_name);
        merged.op.dim = dim;
        merged.op.parity = odd ? Parity::odd : Parity::even;
        ExactMatrix product;
        if (j == i + 1) merged.op.matrix = factors[i].op.matrix;
        else {
            product = multiply_run(factors, i, j, dim);
            if (std::all_of(product.begin(), product.end(), [](auto const& value) { return value.zero(); }))
                return Term{Complex{}, {}, sign};
            needs_scale = !export_matrix(product, merged.op.matrix) || needs_scale;
        }
        // A zero single factor also makes the complete validated term zero.
        if (product.empty() && std::all_of(merged.op.matrix.begin(), merged.op.matrix.end(),
                                         [](auto value) { return value == Complex{}; }))
            return Term{Complex{}, {}, sign};
        out.push_back(std::move(merged));
        products.push_back(std::move(product));
        i = j;
    }

    Term result;
    result.coeff = needs_scale ? export_scaled_products(out, products, t.coeff) : t.coeff;
    for (auto& factor : out) {
        if (matrix_is_identity(factor.op.matrix, factor.op.dim, identity_tolerance)) continue;
        factor.op.qn_deltas = compute_qn_deltas(factor.op.matrix, factor.op.dim,
            factor.basis ? factor.basis->quantum_numbers : std::map<std::string, std::vector<double>>{});
        result.factors.push_back(std::move(factor));
    }
    result.sign = sign;
    return result;
}

TermSum TermSum::located(SiteBasisPtr basis, std::string const& op_name, int site) {
    Term t;
    t.coeff = Complex(1.0);
    Factor f;
    f.site = site;
    f.op = basis->op(op_name);
    f.basis = std::move(basis);
    t.factors.push_back(std::move(f));
    return TermSum({std::move(t)});
}

TermSum TermSum::operator+(TermSum const& other) const {
    std::vector<Term> terms = terms_;
    terms.insert(terms.end(), other.terms_.begin(), other.terms_.end());
    return TermSum(std::move(terms));
}

TermSum& TermSum::operator+=(TermSum const& other) {
    if (&other == this) {
        // Appending a range to itself: reserve first so the source stays valid.
        auto const n = terms_.size();
        terms_.reserve(2 * n);
        for (std::size_t i = 0; i < n; ++i) terms_.push_back(terms_[i]);
        return *this;
    }
    terms_.insert(terms_.end(), other.terms_.begin(), other.terms_.end());
    return *this;
}

TermSum& TermSum::operator+=(TermSum&& other) {
    if (&other == this) return *this += static_cast<TermSum const&>(other);
    if (terms_.empty()) {
        terms_ = std::move(other.terms_);
        return *this;
    }
    terms_.insert(terms_.end(), std::make_move_iterator(other.terms_.begin()),
                  std::make_move_iterator(other.terms_.end()));
    return *this;
}

TermSum& TermSum::operator-=(TermSum const& other) { return *this += -other; }

TermSum TermSum::operator-() const {
    std::vector<Term> terms = terms_;
    for (Term& t : terms) t.sign = -t.sign;
    return TermSum(std::move(terms));
}

TermSum TermSum::operator-(TermSum const& other) const { return *this + (-other); }

TermSum TermSum::operator*(TermSum const& other) const {
    std::vector<Term> terms;
    terms.reserve(terms_.size() * other.terms_.size());
    for (Term const& a : terms_)
        for (Term const& b : other.terms_) {
            Term p;
            p.coeff = a.coeff * b.coeff;
            p.factors = a.factors;
            p.factors.insert(p.factors.end(), b.factors.begin(), b.factors.end());
            p.sign = a.sign * b.sign;
            terms.push_back(std::move(p));
        }
    return TermSum(std::move(terms));
}

TermSum TermSum::operator*(Complex scalar) const {
    if (scalar == Complex(0.0)) return TermSum();
    std::vector<Term> terms = terms_;
    for (Term& t : terms) t.coeff *= scalar;
    return TermSum(std::move(terms));
}

TermSum TermSum::dag() const {
    std::vector<Term> out;
    out.reserve(terms_.size());
    for (Term const& t : terms_) {
        Term d;
        d.coeff = std::conj(t.effective_coeff());
        d.sign = 1;
        d.factors.reserve(t.factors.size());
        for (auto it = t.factors.rbegin(); it != t.factors.rend(); ++it) {
            Factor f;
            f.site = it->site;
            f.op = it->basis ? it->basis->adjoint(it->op) : it->op.dag();
            f.basis = it->basis;
            d.factors.push_back(std::move(f));
        }
        out.push_back(std::move(d));
    }
    return TermSum(std::move(out));
}

TermSum TermSum::canonicalise(bool apply_swap_sign) const {
    return canonicalise(apply_swap_sign, 1e-12);
}

TermSum TermSum::canonicalise(bool apply_swap_sign, double identity_tolerance) const {
    if (!std::isfinite(identity_tolerance) || identity_tolerance < 0)
        throw std::invalid_argument("identity tolerance must be finite and nonnegative");
    // Canonicalise each term, then combine like terms preserving
    // first-appearance order (Python dict insertion order).
    struct Group {
        Term term;
        std::optional<ExactComplex> sum;
        std::vector<Complex> contributions;
    };
    std::vector<Group> merged;
    std::unordered_map<std::string, std::size_t> index;
    for (Term const& raw : terms_) {
        Term t = canonicalise_term(raw, apply_swap_sign, identity_tolerance);
        std::string const key = term_key(t);
        auto it = index.find(key);
        if (it == index.end()) {
            Term entry;
            entry.coeff = t.effective_coeff();
            entry.factors = std::move(t.factors);
            entry.sign = 1;
            index.emplace(key, merged.size());
            merged.push_back({std::move(entry), {}, {}});
        } else {
            auto& group = merged[it->second];
            if (!group.sum) {
                group.sum.emplace(group.term.coeff);
                group.contributions.push_back(group.term.coeff);
            }
            *group.sum += ExactComplex(t.effective_coeff());
            group.contributions.push_back(t.effective_coeff());
        }
    }
    std::vector<Term> out;
    out.reserve(merged.size());
    for (auto& group : merged) {
        if (group.sum && !detail::export_complex(*group.sum, group.term.coeff)) {
            // An overflowing coefficient can still multiply a small finite
            // matrix. Preserve its finite inputs for downstream assembly.
            for (auto coefficient : group.contributions) if (coefficient != Complex{}) {
                out.push_back(group.term);
                out.back().coeff = coefficient;
            }
        } else if (group.term.coeff != Complex{}) out.push_back(std::move(group.term));
    }
    return TermSum(std::move(out));
}

TermSum TermSum::split_factors_by_qn_delta() const {
    std::vector<Term> out;
    for (Term const& term : terms_) {
        // Split components per factor.
        std::vector<std::vector<Factor>> per_factor;
        per_factor.reserve(term.factors.size());
        for (Factor const& f : term.factors) {
            auto const& qns = f.basis ? f.basis->quantum_numbers
                                      : std::map<std::string, std::vector<double>>{};
            if (qns.empty()) { per_factor.push_back({f}); continue; }
            // Inferred deltas use a charge tolerance and cannot certify that
            // all entries belong to one exact component near its boundary.
            std::vector<Factor> choices;
            for (LocalOp& part : f.op.split_by_qn_delta(qns)) {
                Factor comp;
                comp.site = f.site;
                comp.op = std::move(part);
                comp.basis = f.basis;
                choices.push_back(std::move(comp));
            }
            per_factor.push_back(std::move(choices));
        }
        if (per_factor.empty()) { out.push_back(term); continue; }
        // Cartesian product across factors, preserving choice order.
        std::vector<std::vector<Factor>> stack{{}};
        for (auto const& choices : per_factor) {
            std::vector<std::vector<Factor>> next;
            next.reserve(stack.size() * choices.size());
            for (auto const& prefix : stack)
                for (Factor const& choice : choices) {
                    std::vector<Factor> row = prefix;
                    row.push_back(choice);
                    next.push_back(std::move(row));
                }
            stack = std::move(next);
        }
        for (auto& row : stack) {
            Term sub;
            sub.coeff = term.coeff;
            sub.factors = std::move(row);
            sub.sign = term.sign;
            out.push_back(std::move(sub));
        }
    }
    return TermSum(std::move(out));
}

} // namespace models
} // namespace alps
