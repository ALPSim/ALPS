// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// has_sign_problem on kit models: the legacy classifications, per-pair sums
// and inhomogeneous couplings, and the operators outside the rule.
#include <alps/models/builtin.hpp>
#include <alps/models/sign.hpp>
#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace {
namespace lat = alps::lattice;
namespace mdl = alps::models;

mdl::TermSum op(char const* name, int site) { return mdl::TermSum::located(mdl::spin(0.5), name, site); }

mdl::TermSum exchange(int i, int j) {
    return op("Sz", i) * op("Sz", j) + 0.5 * (op("Sp", i) * op("Sm", j) + op("Sm", i) * op("Sp", j));
}

bool ring_with_couplings(std::vector<double> const& couplings) {
    int const n = static_cast<int>(couplings.size());
    mdl::TermSum terms;
    for (int i = 0; i < n; ++i) terms += couplings[std::size_t(i)] * exchange(i, (i + 1) % n);
    return mdl::has_sign_problem(mdl::compile_operator(terms, std::vector<mdl::SiteBasisPtr>(std::size_t(n), mdl::spin(0.5))));
}

TEST(SignProblem, LegacyClassificationsOfSpinModels) {
    EXPECT_FALSE(mdl::has_sign_problem(mdl::heisenberg(lat::chain(8), 1.)));
    EXPECT_TRUE(mdl::has_sign_problem(mdl::heisenberg(lat::chain(7), 1.)));   // odd ring
    EXPECT_FALSE(mdl::has_sign_problem(mdl::heisenberg(lat::chain(7, 1., lat::Boundary::open), 1.)));
    EXPECT_FALSE(mdl::has_sign_problem(mdl::heisenberg(lat::chain(7), -1.)));  // ferromagnet
    EXPECT_FALSE(mdl::has_sign_problem(mdl::heisenberg(lat::square(4, 4), 1.)));
    EXPECT_TRUE(mdl::has_sign_problem(mdl::heisenberg(lat::triangular(3, 3), 1.)));
}

TEST(SignProblem, TransverseFieldsFollowTheSiteSignRule) {
    auto with_field = [](double J) {
        mdl::TermSum terms;
        for (int i = 0; i < 4; ++i) terms += J * exchange(i, (i + 1) % 4) - op("Sx", i);
        return mdl::has_sign_problem(mdl::compile_operator(terms, std::vector<mdl::SiteBasisPtr>(4, mdl::spin(0.5))));
    };
    EXPECT_TRUE(with_field(1.));    // as the legacy example8 antiferromagnet in a transverse field
    EXPECT_FALSE(with_field(-1.));  // and its ferromagnet
}

TEST(SignProblem, SumsEachPairBeforeReadingSigns) {
    // Sx Sx and Sy Sy each have double flips of both signs; their sum has none.
    mdl::TermSum terms;
    for (int i = 0; i < 8; ++i) {
        int const j = (i + 1) % 8;
        terms += op("Sx", i) * op("Sx", j);
        terms += op("Sy", i) * op("Sy", j);
    }
    auto const sites = std::vector<mdl::SiteBasisPtr>(8, mdl::spin(0.5));
    EXPECT_FALSE(mdl::has_sign_problem(mdl::compile_operator(terms, sites)));
    EXPECT_TRUE(mdl::has_sign_problem(mdl::compile_operator(op("Sy", 0) * op("Sy", 1), sites)));
}

TEST(SignProblem, InhomogeneousCouplingsAreCheckedPerPair) {
    EXPECT_FALSE(ring_with_couplings({1, 1, 1, 1, 1, 1, 1, 1}));
    EXPECT_TRUE(ring_with_couplings({1, 1, 1, -1, 1, 1, 1, 1}));    // one ferromagnetic bond frustrates the ring
    EXPECT_FALSE(ring_with_couplings({1, -1, 1, -1, 1, 1, 1, 1}));  // two do not
    EXPECT_FALSE(ring_with_couplings({1, 0, 1, 1, 1, 1, 1}));       // a missing bond opens the odd ring
}

TEST(SignProblem, RejectsOperatorsOutsideTheRule) {
    EXPECT_THROW(mdl::has_sign_problem(mdl::hubbard(lat::chain(4, 1., lat::Boundary::open), 1., 4.)),
                 std::invalid_argument);  // Jordan-Wigner strings
    auto const sites = std::vector<mdl::SiteBasisPtr>(3, mdl::spin(0.5));
    auto const dm = op("Sx", 0) * op("Sy", 1) - op("Sy", 0) * op("Sx", 1);  // Hermitian, imaginary elements
    EXPECT_THROW(mdl::has_sign_problem(mdl::compile_operator(dm, sites)), std::invalid_argument);
    EXPECT_THROW(mdl::has_sign_problem(mdl::compile_operator(op("Sp", 0) * op("Sz", 1) * op("Sm", 2), sites)),
                 std::invalid_argument);
    EXPECT_FALSE(mdl::has_sign_problem(mdl::compile_operator(mdl::TermSum(), sites)));
}

} // namespace
