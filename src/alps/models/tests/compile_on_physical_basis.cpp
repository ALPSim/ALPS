// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// compile_operator against a PhysicalBasis: the same CompiledOperator as
// against the site bases, and a rejection for terms on another local order.
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace {
namespace mdl = alps::models;

void expect_same(mdl::CompiledOperator const& a, mdl::CompiledOperator const& b) {
    ASSERT_EQ(a.site_dims, b.site_dims);
    ASSERT_EQ(a.terms.size(), b.terms.size());
    for (std::size_t t = 0; t < a.terms.size(); ++t) {
        EXPECT_EQ(a.terms[t].coefficient, b.terms[t].coefficient);
        EXPECT_EQ(a.terms[t].jordan_wigner, b.terms[t].jordan_wigner);
        ASSERT_EQ(a.terms[t].factors.size(), b.terms[t].factors.size());
        for (std::size_t f = 0; f < a.terms[t].factors.size(); ++f) {
            EXPECT_EQ(a.terms[t].factors[f].site, b.terms[t].factors[f].site);
            EXPECT_EQ(a.terms[t].factors[f].matrix, b.terms[t].factors[f].matrix);
        }
    }
}

TEST(CompileOnPhysicalBasis, MatchesTheSiteBasesIncludingStrings) {
    auto const model = mdl::hubbard(alps::lattice::chain(5, 1.0, alps::lattice::Boundary::open), 1., 4.);
    auto const physical = mdl::physical_basis(model);
    std::vector<mdl::SiteBasisPtr> sites;
    for (int i = 0; i < model.num_sites(); ++i) sites.push_back(model.site_for(i));
    auto const fermion = model.site_for(0);
    auto located = [&](char const* name, int i) { return mdl::TermSum::located(fermion, name, i); };
    // A string across untouched sites 1-3, a same-site product and a lone odd operator.
    auto const terms = 0.5 * located("c_dag_up", 0) * located("c_up", 4) + located("n_up", 2) * located("n_dn", 2)
                       + located("c_dn", 3);
    expect_same(mdl::compile_operator(terms, physical), mdl::compile_operator(terms, sites, physical.fermion_convention));
    expect_same(mdl::compile_operator(mdl::TermSum(), physical), mdl::compile_operator(mdl::TermSum(), sites));
}

TEST(CompileOnPhysicalBasis, RejectsTermsOnAnotherLocalBasis) {
    auto const model = mdl::heisenberg(alps::lattice::chain(4), 1.);
    auto const physical = mdl::physical_basis(model);
    // Spin-1 terms on a spin-1/2 basis, and a site beyond the basis.
    EXPECT_THROW(mdl::compile_operator(mdl::TermSum::located(mdl::spin(1.), "Sz", 1), physical), std::invalid_argument);
    EXPECT_THROW(mdl::compile_operator(mdl::TermSum::located(mdl::spin(0.5), "Sz", 4), physical), std::invalid_argument);
    // A spin-1/2 site with its states in the other order (as the XML catalog's).
    mdl::SiteBasis reversed = *mdl::spin(0.5);
    reversed.quantum_numbers["Sz"] = {-0.5, 0.5};
    auto const other = mdl::finalize_site(reversed);
    EXPECT_THROW(mdl::compile_operator(mdl::TermSum::located(other, "Sz", 0), physical), std::invalid_argument);
    EXPECT_NO_THROW(mdl::compile_operator(mdl::TermSum::located(mdl::spin(0.5), "Sz", 3), physical));
}

} // namespace
