// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The molecular Hamiltonian factory: term layout, cutoff and input checks.
// Its physics (FCI energies, sign-exact CI vectors) is locked by the
// qdk-chemistry-alps tests against PySCF, through the Python wrapper.
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace mdl = alps::models;

std::vector<std::string> names(mdl::Term const& term) {
    std::vector<std::string> out;
    for (auto const& factor : term.factors) out.push_back(factor.op.name + "@" + std::to_string(factor.site));
    return out;
}

TEST(AbInitio, TermsFollowTheIntegralLoops) {
    int const n = 2;
    std::vector<double> h1 = {-1.25, 0.5, 0.5, 0.0};          // h_11 = 0 is dropped
    std::vector<double> eri(16, 0.);
    eri[0] = 0.75;                                            // (00|00)
    eri[((0 * 2 + 1) * 2 + 1) * 2 + 0] = 1e-13;               // (01|10), below the cutoff
    eri[((1 * 2 + 0) * 2 + 0) * 2 + 1] = 0.25;                // (10|01)
    auto const model = mdl::ab_initio(h1, eri, n);
    EXPECT_EQ(model.num_sites(), 2);
    EXPECT_EQ(model.conserved_qns(), (std::vector<std::string>{"N", "Sz"}));
    auto const terms = model.expand().terms();
    ASSERT_EQ(terms.size(), 3u * 2 + 2u * 4);
    EXPECT_EQ(terms[0].effective_coeff(), mdl::Complex(-1.25));
    EXPECT_EQ(names(terms[0]), (std::vector<std::string>{"c_dag_up@0", "c_up@0"}));
    EXPECT_EQ(names(terms[1]), (std::vector<std::string>{"c_dag_dn@0", "c_dn@0"}));
    EXPECT_EQ(names(terms[2]), (std::vector<std::string>{"c_dag_up@0", "c_up@1"}));
    // 1/2 (pq|rs) c+_{p s1} c+_{r s2} c_{s s2} c_{q s1}, spin pairs (up,up), (up,dn), (dn,up), (dn,dn).
    EXPECT_EQ(terms[6].effective_coeff(), mdl::Complex(0.375));
    EXPECT_EQ(names(terms[7]), (std::vector<std::string>{"c_dag_up@0", "c_dag_dn@0", "c_dn@0", "c_up@0"}));
    EXPECT_EQ(names(terms[10]), (std::vector<std::string>{"c_dag_up@1", "c_dag_up@0", "c_up@1", "c_up@0"}));
    EXPECT_EQ(terms[13].effective_coeff(), mdl::Complex(0.125));
    EXPECT_NO_THROW(mdl::validate_hermiticity(mdl::ab_initio(h1, std::vector<double>(16, 0.), n)));
    EXPECT_EQ(mdl::ab_initio(h1, eri, n, 0.6).expand().terms().size(), 1u * 2 + 1u * 4);  // h_00 and (00|00)
}

TEST(AbInitio, RejectsMismatchedInputs) {
    EXPECT_THROW(mdl::ab_initio({1., 2., 3.}, std::vector<double>(16), 2), std::invalid_argument);
    EXPECT_THROW(mdl::ab_initio(std::vector<double>(4), std::vector<double>(8), 2), std::invalid_argument);
    EXPECT_THROW(mdl::ab_initio(std::vector<double>(4), std::vector<double>(16), 2, -1.), std::invalid_argument);
    EXPECT_THROW(mdl::ab_initio({}, {}, 0), std::invalid_argument);
}

} // namespace
