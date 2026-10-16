// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Equivalent forms of one compiled operator solve alike: factors listed in
// descending site order, scalars moved between coefficients and factor
// matrices, two-site terms expanded into one term per matrix element of the
// first factor (the XML catalog's form), and fermion parity carried by two
// charges (Nup and Ndown) instead of N. Each form is checked against the
// energies of the operator as ALPS::models compiles it.
#include "support.hpp"
#include <gtest/gtest.h>

#include <alps/models/physical_basis.hpp>

#include <algorithm>
#include <complex>
#include <string>
#include <vector>

using namespace classic_test;

namespace {

using Complex = std::complex<double>;

std::vector<double> energies(mdl::CompiledOperator const& op, mdl::PhysicalBasis const& basis, Sector const& sector,
                             std::size_t states, std::size_t count = 1) {
    return cd::solve(op, basis, exact_options(states, sector, count)).energies;
}

void expect_same(std::vector<double> const& actual, std::vector<double> const& expected) {
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t k = 0; k < actual.size(); ++k) EXPECT_NEAR(actual[k], expected[k], 1e-9) << "target " << k;
}

// Every two-site term with its factors in descending site order.
mdl::CompiledOperator reversed(mdl::CompiledOperator op) {
    for (auto& term : op.terms) std::reverse(term.factors.begin(), term.factors.end());
    return op;
}

// The sign of every two-site term moved into its second factor.
mdl::CompiledOperator signs_in_matrices(mdl::CompiledOperator op) {
    for (auto& term : op.terms) {
        if (term.factors.size() != 2) continue;
        term.coefficient = -term.coefficient;
        for (auto& v : term.factors[1].matrix) v = -v;
    }
    return op;
}

// c A x B as sum over (out, in) of (c A[out, in]) |out><in| x B.
mdl::CompiledOperator projector_form(mdl::CompiledOperator const& op) {
    mdl::CompiledOperator out{op.site_dims, {}};
    for (auto const& term : op.terms) {
        if (term.factors.size() != 2) { out.terms.push_back(term); continue; }
        auto const& first = term.factors[0];
        std::size_t const dim = static_cast<std::size_t>(op.site_dims[static_cast<std::size_t>(first.site)]);
        for (std::size_t element = 0; element < first.matrix.size(); ++element) {
            if (first.matrix[element] == Complex{}) continue;
            mdl::LocalTerm piece = term;
            piece.coefficient = term.coefficient * first.matrix[element];
            piece.factors[0].matrix.assign(dim * dim, Complex{});
            piece.factors[0].matrix[element] = 1.;
            out.terms.push_back(std::move(piece));
        }
    }
    return out;
}

} // namespace

TEST(ClassicDmrgOperatorForms, SpinChainFormsAgree) {
    auto const model = mdl::heisenberg(lat::chain(8), 1.0, 0.0, 1.0); // periodic spin-1 ring
    auto const op = mdl::compile_hamiltonian(model).op;
    auto const basis = mdl::physical_basis(model);
    Sector const sector{{"Sz", 0.}};
    auto const expected = ed_energies(model, sector, 2);
    expect_same(energies(op, basis, sector, 200, 2), expected);
    expect_same(energies(reversed(op), basis, sector, 200, 2), expected);
    expect_same(energies(signs_in_matrices(op), basis, sector, 200, 2), expected);
    expect_same(energies(projector_form(op), basis, sector, 200, 2), expected);
}

TEST(ClassicDmrgOperatorForms, FermionFormsAgree) {
    auto const model = mdl::hubbard(lat::chain(6), 1.0, 4.0); // periodic ring: wrap strings
    auto const op = mdl::compile_hamiltonian(model).op;
    auto const basis = mdl::physical_basis(model);
    Sector const sector{{"N", 4.}, {"Sz", 0.}};
    auto const expected = ed_energies(model, sector, 2);
    expect_same(energies(op, basis, sector, 256, 2), expected);
    expect_same(energies(reversed(op), basis, sector, 256, 2), expected);
    expect_same(energies(signs_in_matrices(op), basis, sector, 256, 2), expected);
    expect_same(energies(projector_form(reversed(op)), basis, sector, 256, 2), expected);
}

TEST(ClassicDmrgOperatorForms, TwoChargesCarryTheFermionParity) {
    auto const model = mdl::hubbard(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 4.0);
    auto const op = mdl::compile_hamiltonian(model).op;
    // Nup and Ndown in place of N and Sz: neither alone has the state parity.
    auto basis = mdl::physical_basis(model);
    for (auto& site : basis.sites) {
        auto const n = site.quantum_numbers.at("N"), sz = site.quantum_numbers.at("Sz");
        std::vector<double> up, down;
        for (std::size_t s = 0; s < n.size(); ++s) {
            up.push_back(n[s] / 2 + sz[s]);
            down.push_back(n[s] / 2 - sz[s]);
        }
        site.quantum_numbers = {{"Nup", up}, {"Ndown", down}};
    }
    basis.conserved_quantum_numbers = {"Ndown", "Nup"};
    expect_same(energies(op, basis, {{"Nup", 3.}, {"Ndown", 2.}}, 256, 2),
                ed_energies(model, {{"N", 5.}, {"Sz", 0.5}}, 2));
    expect_same(energies(op, basis, {{"Nup", 3.}, {"Ndown", 3.}}, 256),
                ed_energies(model, {{"N", 6.}, {"Sz", 0.}}));
    // With no charge conserved the parity still needs both.
    basis.conserved_quantum_numbers.clear();
    EXPECT_NO_THROW(energies(op, basis, {}, 64));
    for (auto& site : basis.sites) site.quantum_numbers.erase("Ndown");
    EXPECT_THROW(energies(op, basis, {}, 64), std::invalid_argument);
}
