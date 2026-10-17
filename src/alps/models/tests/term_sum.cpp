// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// TermSum in-place arithmetic and the expansion order built on it.
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <complex>
#include <utility>
#include <vector>

namespace {
namespace mdl = alps::models;
using C = mdl::Complex;

void expect_same_terms(mdl::TermSum const& actual, mdl::TermSum const& expected) {
    ASSERT_EQ(actual.terms().size(), expected.terms().size());
    for (std::size_t t = 0; t < expected.terms().size(); ++t) {
        auto const& a = actual.terms()[t];
        auto const& e = expected.terms()[t];
        EXPECT_EQ(a.coeff, e.coeff) << "term " << t;
        EXPECT_EQ(a.sign, e.sign) << "term " << t;
        ASSERT_EQ(a.factors.size(), e.factors.size()) << "term " << t;
        for (std::size_t f = 0; f < e.factors.size(); ++f) {
            EXPECT_EQ(a.factors[f].site, e.factors[f].site);
            EXPECT_EQ(a.factors[f].op.name, e.factors[f].op.name);
            EXPECT_EQ(a.factors[f].op.matrix, e.factors[f].op.matrix);
            EXPECT_EQ(a.factors[f].basis, e.factors[f].basis);
        }
    }
}

mdl::TermSum op(mdl::SiteBasisPtr const& site, char const* name, int i) {
    return mdl::TermSum::located(site, name, i);
}

TEST(TermSum, InPlaceSumAppendsLikeTheBinarySum) {
    auto const spin = mdl::spin(0.5);
    auto const a = 0.5 * op(spin, "Sp", 0) * op(spin, "Sm", 1) + C(0., 2.) * op(spin, "Sz", 2);
    auto const b = -1.5 * op(spin, "Sz", 0) * op(spin, "Sz", 3) - op(spin, "Sx", 1);

    auto sum = a;
    sum += b;
    expect_same_terms(sum, a + b);

    auto difference = a;
    difference -= b;
    expect_same_terms(difference, a - b);

    auto moved = a;
    auto source = b;
    moved += std::move(source);
    expect_same_terms(moved, a + b);

    mdl::TermSum empty;
    empty += mdl::TermSum(b);
    expect_same_terms(empty, b);
}

TEST(TermSum, InPlaceSumWithItselfDoublesTheTerms) {
    auto const spin = mdl::spin(0.5);
    auto const a = op(spin, "Sz", 0) + 2. * op(spin, "Sp", 1) * op(spin, "Sm", 2);
    auto doubled = a;
    doubled += doubled;
    expect_same_terms(doubled, a + a);
    auto moved = a;
    moved += std::move(moved);
    expect_same_terms(moved, a + a);
    auto none = a;
    none -= none;
    expect_same_terms(none, a - a);
}

// expand() keeps its template-major, bond/site-minor order and places each
// Hermitian conjugate right after the contribution it conjugates.
TEST(TermSum, ExpansionOrderMatchesTheBinarySumReference) {
    auto const lattice = alps::lattice::chain(4);
    auto const spin = mdl::spin(0.5);
    mdl::Hamiltonian model(lattice, spin, mdl::FermionConvention::strict, {"Sz"});
    model.add_bond(0.5, [&](int i, int j) { return op(spin, "Sp", i) * op(spin, "Sm", j); }, std::nullopt, true);
    model.add_site(C(0.3), [&](int i) { return op(spin, "Sz", i); });
    auto const global = 0.25 * op(spin, "Sp", 0) * op(spin, "Sm", 2);
    model.add(global, true);
    model.add(op(spin, "Sz", 1) * op(spin, "Sz", 3));

    mdl::TermSum reference;
    for (auto const& bond : lattice.bonds()) {
        auto const term = op(spin, "Sp", bond.source) * op(spin, "Sm", bond.target) * C(0.5);
        reference = reference + (term + term.dag());
    }
    for (auto const& site : lattice.sites()) reference = reference + op(spin, "Sz", site.index) * C(0.3);
    reference = reference + (global + global.dag());
    reference = reference + op(spin, "Sz", 1) * op(spin, "Sz", 3);
    expect_same_terms(model.expand(), reference);
}

} // namespace
