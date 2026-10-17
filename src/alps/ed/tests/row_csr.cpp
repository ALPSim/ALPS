// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// build_canonical_csr (row-wise, from the adjoint) against the COO path it
// replaces, canonical_csr(build_coo_from_local): the same real/complex choice
// and bit-identical offsets, indices and values.
#include <alps/ed/product_basis.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
namespace ed = alps::ed;
namespace mdl = alps::models;
using C = std::complex<double>;

std::vector<std::int64_t> sector_states(std::vector<int> const& dims, std::vector<std::vector<double>> const& charge,
                                        double target) {
    std::vector<ed::BasisConstraint> constraints;
    if (!charge.empty()) {
        auto& c = constraints.emplace_back();
        c.target = target;
        c.values = charge;
    }
    auto const basis = ed::enumerate_product_basis(dims, constraints);
    return {basis.states.begin(), basis.states.end()};
}

void expect_same_csr(mdl::CompiledOperator const& op, std::vector<std::int64_t> const& states,
                     std::vector<std::vector<double>> const& parity, bool expect_complex) {
    auto const count = op.site_dims.empty() ? 0 : states.size() / op.site_dims.size();
    auto const coo = ed::build_coo_from_local(op, states, count, parity);
    ASSERT_EQ(coo.is_complex, expect_complex);
    auto const rows = ed::build_canonical_csr(op, states, count, parity);
    ASSERT_EQ(rows.index(), coo.is_complex ? 1u : 0u);
    auto const check = [&](auto const& actual, auto const& expected) {
        EXPECT_EQ(actual.rows, expected.rows);
        EXPECT_EQ(actual.columns, expected.columns);
        EXPECT_EQ(actual.offsets, expected.offsets);
        EXPECT_EQ(actual.indices, expected.indices);
        ASSERT_EQ(actual.values.size(), expected.values.size());
        for (std::size_t i = 0; i < expected.values.size(); ++i)
            ASSERT_EQ(actual.values[i], expected.values[i]) << "entry " << i;  // bit for bit, not to a tolerance
    };
    if (coo.is_complex) check(std::get<1>(rows), ed::canonical_csr<C>(coo, count, count));
    else check(std::get<0>(rows), ed::canonical_csr<double>(coo, count, count));
    // Forced complex output equals the complex canonical CSR either way.
    check(std::get<1>(ed::build_canonical_csr(op, states, count, parity, true)), ed::canonical_csr<C>(coo, count, count));
}

std::vector<std::vector<double>> parities(std::vector<mdl::SiteBasisPtr> const& sites) {
    std::vector<std::vector<double>> out;
    for (auto const& site : sites) out.push_back(site->parity_operator);
    return out;
}

TEST(RowCsr, RealSpinSector) {
    auto const model = mdl::heisenberg(alps::lattice::chain(8), 1., 0.3);
    auto const op = mdl::compile_hamiltonian(model).op;
    std::vector<std::vector<double>> sz(8, {0.5, -0.5});
    expect_same_csr(op, sector_states(op.site_dims, sz, 0.), {}, false);
    expect_same_csr(op, sector_states(op.site_dims, {}, 0.), {}, false);
}

TEST(RowCsr, ComplexTermsAndNonHermitianOperators) {
    auto const spin = mdl::spin(0.5);
    std::vector<mdl::SiteBasisPtr> const sites(5, spin);
    auto located = [&](char const* name, int i) { return mdl::TermSum::located(spin, name, i); };
    mdl::TermSum terms;
    for (int i = 0; i + 1 < 5; ++i) {
        terms += located("Sx", i) * located("Sy", i + 1) - located("Sy", i) * located("Sx", i + 1);  // DM, complex
        terms += 0.7 * located("Sz", i) * located("Sz", i + 1);
    }
    terms += C(0., 0.4) * located("Sy", 2);
    auto const op = mdl::compile_operator(terms, sites);
    auto const all = sector_states(op.site_dims, {}, 0.);
    expect_same_csr(op, all, {}, true);
    // Non-Hermitian and charge-changing: a lone hop and a lone raising operator.
    expect_same_csr(mdl::compile_operator(1.5 * located("Sp", 0) * located("Sm", 3) + located("Sp", 4), sites), all, {}, false);
}

TEST(RowCsr, ComplexEntriesThatNeverContributeAreReal) {
    auto const spin = mdl::spin(0.5);
    std::vector<mdl::SiteBasisPtr> const sites(4, spin);
    auto located = [&](char const* name, int i) { return mdl::TermSum::located(spin, name, i); };
    // i Sp leaves the Sz=0 sector: no contribution, so the COO path calls it real.
    auto const op = mdl::compile_operator(located("Sz", 0) * located("Sz", 1) + C(0., 1.) * located("Sp", 2), sites);
    std::vector<std::vector<double>> sz(4, {0.5, -0.5});
    expect_same_csr(op, sector_states(op.site_dims, sz, 0.), {}, false);
}

TEST(RowCsr, FermionStringsOnFactorAndSpectatorSites) {
    auto const site = mdl::spinful_fermion();
    int const n = 5;
    std::vector<mdl::SiteBasisPtr> const sites(n, site);
    auto c = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-1., 1.);
    mdl::TermSum terms;
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            auto const t = u(rng);
            terms += t * (c("c_dag_up", p) * c("c_up", q) + c("c_dag_dn", p) * c("c_dn", q));
        }
    for (int k = 0; k < 12; ++k) {  // four-operator terms on random orbitals, both spins
        int const p = rng() % n, q = rng() % n, r = rng() % n, s = rng() % n;
        terms += u(rng) * c("c_dag_up", p) * c("c_dag_dn", r) * c("c_dn", s) * c("c_up", q);
    }
    auto const op = mdl::compile_operator(terms, sites);
    ASSERT_TRUE(std::any_of(op.terms.begin(), op.terms.end(), [](auto const& term) {
        return std::any_of(term.factors.begin(), term.factors.end(), [&](auto const& f) {
            return std::count(term.jordan_wigner.begin(), term.jordan_wigner.end(), f.site) > 0;
        });
    })) << "the case must place a string on a factor site";
    std::vector<std::vector<double>> number(n, {0., 1., 1., 2.});
    expect_same_csr(op, sector_states(op.site_dims, number, 5.), parities(sites), false);
    expect_same_csr(op, sector_states(op.site_dims, number, 4.), {}, false);  // F = 1 without parity
}

TEST(RowCsr, DuplicatesSumInTermOrderAndExactCancellationsVanish) {
    auto const sz = mdl::spin(0.5)->op("Sz").matrix;
    auto const sp = mdl::spin(0.5)->op("Sp").matrix;
    auto const sm = mdl::spin(0.5)->op("Sm").matrix;
    auto term = [](C coefficient, std::vector<mdl::LocalFactor> factors) {
        mdl::LocalTerm t;
        t.coefficient = coefficient;
        t.factors = std::move(factors);
        return t;
    };
    mdl::CompiledOperator op;
    op.site_dims = {2, 2, 2};
    // Uncombined like terms: an exact cancellation and sums whose rounding depends on order.
    op.terms = {term(0.5, {{0, sz}, {1, sz}}), term(1e16, {{0, sp}, {1, sm}}), term(-0.5, {{0, sz}, {1, sz}}),
                term(1., {{0, sp}, {1, sm}}), term(-1e16, {{0, sp}, {1, sm}}), term(0.1, {{2, sz}}),
                term(0.2, {{2, sz}}), term(0.3, {{2, sz}})};
    expect_same_csr(op, sector_states(op.site_dims, {}, 0.), {}, false);
}

TEST(RowCsr, UnsortedAndUnencodableBases) {
    auto const spin = mdl::spin(0.5);
    std::vector<mdl::SiteBasisPtr> sites(6, spin);
    auto located = [&](char const* name, int i) { return mdl::TermSum::located(spin, name, i); };
    mdl::TermSum terms;
    for (int i = 0; i + 1 < 6; ++i) terms += located("Sp", i) * located("Sm", i + 1) + located("Sm", i) * located("Sp", i + 1);
    auto const op = mdl::compile_operator(terms, sites);
    std::vector<std::vector<double>> sz(6, {0.5, -0.5});
    auto states = sector_states(op.site_dims, sz, 0.);
    // Reverse the rows: the lookup table path, rows kept in the caller's order.
    std::vector<std::int64_t> reversed;
    for (std::size_t row = states.size() / 6; row-- > 0;) reversed.insert(reversed.end(), states.begin() + row * 6, states.begin() + row * 6 + 6);
    expect_same_csr(op, reversed, {}, false);

    // 70 two-level sites do not fit uint64 codes: the explicit-row index.
    int const n = 70;
    sites.assign(n, spin);
    mdl::TermSum chain;
    for (int i = 0; i + 1 < n; ++i) chain += located("Sp", i) * located("Sm", i + 1) + located("Sm", i) * located("Sp", i + 1);
    auto const wide = mdl::compile_operator(chain, sites);
    std::vector<std::int64_t> one_flip;
    for (int k = 0; k < n; ++k)
        for (int i = 0; i < n; ++i) one_flip.push_back(i == k ? 0 : 1);
    expect_same_csr(wide, one_flip, {}, false);
}

TEST(RowCsr, RejectsWhatTheCooPathRejects) {
    auto const spin = mdl::spin(0.5);
    std::vector<mdl::SiteBasisPtr> const sites(2, spin);
    auto const op = mdl::compile_operator(mdl::TermSum::located(spin, "Sz", 0), sites);
    EXPECT_THROW(ed::build_canonical_csr(op, {0, 0, 0, 0}, 2), std::invalid_argument);   // repeated state
    EXPECT_THROW(ed::build_canonical_csr(op, {0, 2}, 1), std::invalid_argument);         // local state out of range
    EXPECT_THROW(ed::build_canonical_csr(op, {0, 1, 1}, 2), std::invalid_argument);      // shape mismatch
    EXPECT_THROW(ed::build_canonical_csr(op, {0, 1}, 1, {{1., -1.}}), std::invalid_argument);  // parity site count
    auto bad = op;
    bad.terms[0].factors[0].site = 5;
    EXPECT_THROW(ed::build_canonical_csr(bad, {0, 1}, 1), std::invalid_argument);
    bad = op;
    bad.terms[0].coefficient = C(std::numeric_limits<double>::infinity());
    EXPECT_THROW(ed::build_canonical_csr(bad, {0, 1}, 1), std::invalid_argument);
}

} // namespace
