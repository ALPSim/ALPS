// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Catalog spellings of the built-in operators, and adjoints named by the
// site that defines them.
#include <alps/models/operator_names.hpp>
#include <alps/models/site.hpp>
#include <alps/models/terms.hpp>
#include <gtest/gtest.h>

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace mdl = alps::models;

mdl::Matrix dagger(mdl::LocalOp const& op) {
    mdl::Matrix out(op.matrix.size());
    for (int r = 0; r < op.dim; ++r)
        for (int c = 0; c < op.dim; ++c) out[static_cast<std::size_t>(c) * op.dim + r] = std::conj(op(r, c));
    return out;
}

std::vector<mdl::SiteBasisPtr> built_in_sites() {
    return {mdl::spin(0.5), mdl::spin(1.5), mdl::boson(3), mdl::hardcore_boson(), mdl::spinless_fermion(),
            mdl::spinful_fermion()};
}

TEST(OperatorNames, CatalogSpellingsReachTheBuiltInOperators) {
    auto const spin = mdl::spin(0.5);
    EXPECT_EQ(spin->op("Splus").name, "Sp");
    EXPECT_EQ(spin->op("Sminus").matrix, spin->op("Sm").matrix);
    EXPECT_TRUE(spin->has_op("Splus"));
    EXPECT_FALSE(spin->has_op("cdag_up"));
    EXPECT_THROW(spin->op("cdag_up"), std::out_of_range);
    auto const fermion = mdl::spinful_fermion();
    EXPECT_EQ(fermion->op("cdag_up").name, "c_dag_up");
    EXPECT_EQ(fermion->op("cdag_down").name, "c_dag_dn");
    EXPECT_EQ(fermion->op("c_down").name, "c_dn");
    EXPECT_EQ(fermion->op("n_down").name, "n_dn");
    EXPECT_EQ(fermion->op("double_occupancy").name, "n_updn");
    EXPECT_EQ(fermion->op("Splus").name, "Sp");
    EXPECT_EQ(mdl::spinless_fermion()->op("cdag").name, "c_dag");
    EXPECT_EQ(mdl::boson(2)->op("bdag").name, "b_dag");  // the first kit name the site defines
    EXPECT_EQ(mdl::boson(2)->op("b").name, "b");          // a direct name wins over an alias

    // Every kit name in the table is a built-in operator somewhere.
    for (auto const& [kit, catalog] : mdl::catalog_operator_names) {
        bool defined = false;
        for (auto const& site : built_in_sites()) defined = defined || site->operators.count(std::string(kit));
        EXPECT_TRUE(defined) << kit;
    }
    // located() takes either spelling and stores the kit operator.
    auto const term = mdl::TermSum::located(fermion, "cdag_up", 2);
    EXPECT_EQ(term.terms()[0].factors[0].op.name, "c_dag_up");
}

TEST(OperatorNames, AdjointsCarryTheSitesOwnNames) {
    auto const spin = mdl::spin(0.5);
    EXPECT_EQ(spin->adjoint(spin->op("Sp")).name, "Sm");
    EXPECT_EQ(spin->adjoint(spin->op("Sm")).name, "Sp");
    for (char const* hermitian : {"Sz", "Sx", "Sy", "S2", "I"}) EXPECT_EQ(spin->adjoint(spin->op(hermitian)).name, hermitian);
    auto const fermion = mdl::spinful_fermion();
    EXPECT_EQ(fermion->adjoint(fermion->op("c_up")).name, "c_dag_up");
    EXPECT_EQ(fermion->adjoint(fermion->op("c_dag_dn")).name, "c_dn");
    EXPECT_EQ(fermion->adjoint(fermion->op("n_updn")).name, "n_updn");
    auto const boson = mdl::boson(3);
    EXPECT_EQ(boson->adjoint(boson->op("a")).name, "a_dag");
    EXPECT_EQ(boson->adjoint(boson->op("b_dag")).name, "b");

    // Whatever name comes back, the matrix is the adjoint, and a named one is the site's own.
    for (auto const& site : built_in_sites())
        for (auto const& [name, op] : site->operators) {
            auto const adjoint = site->adjoint(op);
            EXPECT_EQ(adjoint.matrix, dagger(op)) << name;
            EXPECT_EQ(adjoint.parity, op.parity) << name;
            ASSERT_TRUE(site->operators.count(adjoint.name)) << name << " -> " << adjoint.name;
            EXPECT_EQ(site->operators.at(adjoint.name).matrix, adjoint.matrix) << name;
        }

    // An operator the site does not define keeps the suffix rule.
    auto custom = spin->op("Sp");
    custom.matrix[1] = mdl::Complex(2.);
    custom.name = "twice_Sp";
    EXPECT_EQ(spin->adjoint(custom).name, "twice_Sp_dag");
}

TEST(OperatorNames, ConjugatedTermSumsNameTheirFactors) {
    auto const spin = mdl::spin(0.5);
    auto const hop = mdl::TermSum::located(spin, "Sp", 0) * mdl::TermSum::located(spin, "Sm", 3);
    auto const conjugate = hop.dag();
    ASSERT_EQ(conjugate.terms().size(), 1u);
    auto const& factors = conjugate.terms()[0].factors;
    ASSERT_EQ(factors.size(), 2u);
    EXPECT_EQ(factors[0].op.name, "Sp");
    EXPECT_EQ(factors[0].site, 3);
    EXPECT_EQ(factors[1].op.name, "Sm");
    EXPECT_EQ(factors[1].site, 0);
    auto const fermion = mdl::spinful_fermion();
    auto const pair = (mdl::TermSum::located(fermion, "c_dag_up", 1) * mdl::TermSum::located(fermion, "c_dn", 2)).dag();
    EXPECT_EQ(pair.terms()[0].factors[0].op.name, "c_dag_dn");
    EXPECT_EQ(pair.terms()[0].factors[1].op.name, "c_up");
}

} // namespace
