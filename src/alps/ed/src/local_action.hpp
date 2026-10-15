// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/assembly/local_action.hpp.
// Adaptation: <alps/ed/...> include names.
#pragma once
#include <alps/models/compiled.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
namespace alps::ed::detail {
// Precompute nonzero outputs for each local input column once per term.
// Both language-facing assembly paths use this representation and scatter.
template <typename Scalar> struct PreparedTerm {
    struct Action {
        int site;
        std::vector<std::vector<std::pair<int, Scalar>>> columns;
    };
    Scalar coeff;
    std::vector<Action> actions;
    std::vector<int> jw;
    bool single_output = true;
};

template <typename Scalar, typename Convert>
std::vector<PreparedTerm<Scalar>> prepare_terms(models::CompiledOperator const &input, Convert convert,
                                                bool require_unique_sites) {
    auto const n = input.site_dims.size();
    auto const count = input.terms.size();
    std::vector<PreparedTerm<Scalar>> terms;
    terms.reserve(count);
    for (auto const &source : input.terms) {
        PreparedTerm<Scalar> term{convert(source.coefficient, "coefficient"), {}, source.jordan_wigner, true};
        std::vector<bool> seen(n, false);
        for (auto const &factor : source.factors) {
            int const site = factor.site;
            if (site < 0 || static_cast<std::size_t>(site) >= n)
                throw std::invalid_argument("sparse assembly: factor site out of range");
            if (require_unique_sites && seen[static_cast<std::size_t>(site)])
                throw std::invalid_argument("sparse assembly: compose repeated factor sites first");
            seen[static_cast<std::size_t>(site)] = true;
            auto const d = static_cast<std::size_t>(input.site_dims[static_cast<std::size_t>(site)]);
            auto const &matrix = factor.matrix;
            if (d > std::numeric_limits<std::size_t>::max() / d || matrix.size() != d * d)
                throw std::invalid_argument("sparse assembly: operator size != site dim^2");
            typename PreparedTerm<Scalar>::Action action{site, {}};
            action.columns.resize(d);
            for (std::size_t row = 0; row < d; ++row)
                for (std::size_t col = 0; col < d; ++col) {
                    auto const value = convert(matrix[row * d + col], "operator entry");
                    if (value != Scalar{})
                        action.columns[col].emplace_back(static_cast<int>(row), value);
                }
            for (auto const &column : action.columns)
                if (column.size() > 1)
                    term.single_output = false;
            term.actions.push_back(std::move(action));
        }
        for (int site : term.jw)
            if (site < 0 || static_cast<std::size_t>(site) >= n)
                throw std::invalid_argument("sparse assembly: JW site out of range");
        terms.push_back(std::move(term));
    }
    return terms;
}

template <typename Scalar, typename Basis, typename Emit> struct TermScatter {
    PreparedTerm<Scalar> const &term;
    Basis const &basis;
    std::size_t input;
    Emit &emit;

    void descend(std::size_t k, Scalar amplitude, typename Basis::State &target) const {
        if (k == term.actions.size()) {
            emit(target, amplitude);
            return;
        }
        auto const &action = term.actions[k];
        int const before = basis.digit(input, action.site);
        for (auto const &[after, value] : action.columns[static_cast<std::size_t>(before)]) {
            basis.shift(target, action.site, before, after);
            descend(k + 1, amplitude * value, target);
            basis.shift(target, action.site, after, before);
        }
    }
};

template <typename Scalar, typename Basis, typename Emit>
void apply_term_state(PreparedTerm<Scalar> const &term, Basis const &basis,
                      std::vector<std::vector<double> const *> const &parity, std::size_t input, Emit emit) {
    double sign = 1.;
    for (int site : term.jw)
        if (parity[static_cast<std::size_t>(site)])
            sign *= (*parity[static_cast<std::size_t>(site)])[static_cast<std::size_t>(basis.digit(input, site))];
    if (sign == 0.)
        return;
    Scalar const amplitude = term.coeff * sign;
    if (term.actions.empty()) {
        emit(basis.state(input), amplitude);
        return;
    }
    auto target = basis.state(input);
    if (term.single_output) {
        // Diagonal, raising/lowering and permutation factors cannot branch.
        // Preserve the same multiplication order without recursive scatter.
        Scalar value = amplitude;
        for (auto const &action : term.actions) {
            int const before = basis.digit(input, action.site);
            auto const &column = action.columns[static_cast<std::size_t>(before)];
            if (column.empty())
                return;
            auto const &[after, entry] = column.front();
            basis.shift(target, action.site, before, after);
            value *= entry;
        }
        emit(target, value);
    } else {
        TermScatter<Scalar, Basis, Emit>{term, basis, input, emit}.descend(0, amplitude, target);
    }
}

template <typename Scalar, typename Basis, typename Emit>
void apply_term(PreparedTerm<Scalar> const &term, Basis const &basis,
                std::vector<std::vector<double> const *> const &parity, std::size_t input, Emit emit) {
    apply_term_state(term, basis, parity, input, [&](auto const &target, Scalar amplitude) {
        auto const output =
            basis.unchanged(target, input) ? static_cast<std::int64_t>(input) : basis.index(target);
        if (output >= 0)
            emit(output, amplitude);
    });
}

} // namespace alps::ed::detail
