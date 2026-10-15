// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/model/tests/models/catalog_compiler.cpp. Adaptation: the shipped XML
// catalogs replace alpsx's TOML fixtures ("boson Hubbard" with Nmax=1 and
// "spinless fermions" stand in for its "boson" and "complex fermion" models);
// the charge accessors are new.
// This translation unit deliberately includes no ED or matrix-storage header.
#include <alps/model/catalog.hpp>
#include <gtest/gtest.h>

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
alps::Parameters parameters(int sites) {
    alps::Parameters p;
    p["MODEL"] = "boson Hubbard";
    p["LATTICE"] = "chain lattice";
    p["L"] = sites;
    p["Nmax"] = 1;
    p["N_total"] = 1;
    p["t"] = .7;
    return p;
}
} // namespace

TEST(CatalogCompiler, LocalCompilationDoesNotEnumerateTheProductSpace) {
    // The unconstrained local space is 2^1000. Even the two-particle sector
    // would have 499,500 states; the model adapter builds only local metadata.
    auto p = parameters(1000);
    p["N_total"] = 2;
    alps::models::CatalogModel<> model(p);
    EXPECT_EQ(model.basis().size(), 1000u);
    auto const op = model.compile_hamiltonian();
    EXPECT_EQ(op.site_dims, std::vector<int>(1000, 2));
    EXPECT_FALSE(op.terms.empty());
    EXPECT_LT(op.terms.size(), 10000u);
    EXPECT_EQ(model.local_parity(), std::vector<std::vector<double>>(1000, {1., 1.}));
    EXPECT_EQ(model.conserved_quantum_numbers(), std::vector<std::string>{"N"});
    EXPECT_EQ(model.sector(), (std::map<std::string, double>{{"N", 2.}}));
    EXPECT_EQ(model.local_charge("N"), std::vector<std::vector<double>>(1000, {0., 1.}));
}

TEST(CatalogCompiler, InputsAndCompiledOperatorsOwnTheirData) {
    auto p = parameters(3);
    alps::models::CatalogModel<> model(p);
    auto original = model.compile_hamiltonian();
    auto const original_value = original;
    p["N_total"] = 2;
    EXPECT_EQ(alps::to_double(model.basis().get_basis().constraints()[0].second), 1.);
    model.set_parameters(p);
    EXPECT_EQ(alps::to_double(model.basis().get_basis().constraints()[0].second), 2.);
    // A charge restriction affects enumeration, not the local operator.
    EXPECT_EQ(model.compile_hamiltonian(), original);
    original.terms.front().coefficient = 17.;
    EXPECT_EQ(model.compile_hamiltonian(), original_value);
    EXPECT_NE(model.compile_hamiltonian(), original);
    const auto detached = [] {
        alps::models::CatalogModel<> temporary(parameters(3));
        return temporary.compile_hamiltonian();
    }();
    EXPECT_EQ(detached, original_value);
}

TEST(CatalogCompiler, ObservableOrientationAndParityAreExplicit) {
    auto p = parameters(3);
    p["MODEL"] = "spinless fermions";
    alps::models::CatalogModel<> model(p);
    EXPECT_EQ(model.local_parity(), std::vector<std::vector<double>>(3, {1., -1.}));
    // Catalog site tensors reject odd fermionic operators; this adapter
    // retains that restriction rather than inventing a parity string.
    EXPECT_THROW(model.compile_operator(alps::SiteOperator("c(i)", "i"),
                                       decltype(model)::site_descriptor(1)), std::runtime_error);
    alps::models::CatalogModel<> boson(parameters(3));
    auto const op = boson.compile_operator(alps::SiteOperator("b(i)", "i"),
                                         decltype(boson)::site_descriptor(1));
    ASSERT_EQ(op.terms.size(), 1u);
    ASSERT_EQ(op.terms[0].factors.size(), 1u);
    EXPECT_EQ(op.terms[0].factors[0].site, 1);
    // OUT-row/IN-column: b takes the occupied state (column 1) to the vacuum (row 0).
    EXPECT_EQ(op.terms[0].factors[0].matrix,
              (std::vector<std::complex<double>>{0., 1., 0., 0.}));
}

TEST(CatalogCompiler, ChargesComeFromTheBasisConstraints) {
    auto p = parameters(3);
    p["MODEL"] = "fermion Hubbard";
    p.erase("N_total");
    p["Nup_total"] = 1;
    p["CONSERVED_QUANTUMNUMBERS"] = "Ndown, Nup";
    alps::models::CatalogModel<> model(p);
    EXPECT_EQ(model.conserved_quantum_numbers(), (std::vector<std::string>{"Ndown", "Nup"}));
    EXPECT_EQ(model.sector(), (std::map<std::string, double>{{"Nup", 1.}}));
    // Local states in catalog order: (Nup, Ndown) = (0,0), (0,1), (1,0), (1,1).
    EXPECT_EQ(model.local_charge("Nup"), std::vector<std::vector<double>>(3, {0., 0., 1., 1.}));
    EXPECT_EQ(model.local_parity(), std::vector<std::vector<double>>(3, {1., -1., -1., 1.}));
    EXPECT_THROW(model.local_charge("Sz"), std::invalid_argument);
    p["CONSERVED_QUANTUMNUMBERS"] = "Sz";
    alps::models::CatalogModel<> undeclared(p);
    EXPECT_THROW(undeclared.conserved_quantum_numbers(), std::invalid_argument);
}
