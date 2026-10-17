// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Rejected inputs of the classic DMRG solver: option ranges, sector
// declarations and feasibility, term shapes the engine cannot represent,
// non-Hermitian or non-conserving operators and recursive invocation.
#include "support.hpp"
#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

using namespace classic_test;

namespace {
mdl::Hamiltonian chain(int length = 6) { return mdl::heisenberg(lat::chain(length, 1.0, lat::Boundary::open), 1.0); }
} // namespace

TEST(ClassicDmrgValidation, SweepSchedulesAreChecked) {
    EXPECT_THROW(cd::fixed_sweeps(0, 10), std::invalid_argument);
    EXPECT_THROW(cd::fixed_sweeps(2, 0), std::invalid_argument);
    EXPECT_THROW(cd::ramped_sweeps(4, 7), std::invalid_argument);
    auto const ramp = cd::ramped_sweeps(2, 8);
    ASSERT_EQ(ramp.size(), 2u);
    EXPECT_EQ(ramp[0], (cd::SweepStates{2, 4}));
    EXPECT_EQ(ramp[1], (cd::SweepStates{6, 8}));
    cd::Options options = exact_options(8, {{"Sz", 0.}});
    options.sweeps.clear();
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
    options = exact_options(8, {{"Sz", 0.}});
    options.sweeps[1].left_to_right = 0;
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
}

TEST(ClassicDmrgValidation, ControlRangesAreChecked) {
    auto base = exact_options(8, {{"Sz", 0.}});
    auto options = base;
    options.warmup_states = 0;
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
    options = base;
    options.num_eigenvalues = 0;
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
    options = base;
    options.truncation_error = 1.5;
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
    options = base;
    options.truncation_max_states = 10;
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
    options = base;
    options.lanczos_tolerance = -1.;
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
    options = base;
    options.scratch_parent = "/nonexistent/classic-dmrg-scratch";
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument);
}

TEST(ClassicDmrgValidation, SectorsMustBeDeclaredRepresentableAndNonempty) {
    EXPECT_THROW(cd::solve(chain(), exact_options(8, {{"N", 0.}})), std::invalid_argument);
    EXPECT_THROW(cd::solve(chain(), exact_options(8, {{"Sz", 0.3}})), std::invalid_argument);
    EXPECT_THROW(cd::solve(chain(), exact_options(8, {{"Sz", 5.}})), std::invalid_argument);
    EXPECT_THROW(cd::solve(chain(), exact_options(8, {{"Sz", 0.5}})), std::invalid_argument); // parity mismatch
    EXPECT_THROW(cd::solve(chain(4), exact_options(8, {{"Sz", 0.}}, 100)), std::invalid_argument);
    EXPECT_THROW(cd::solve(chain(4), exact_options(8, {}, 100)), std::invalid_argument);
    EXPECT_THROW(cd::solve(chain(3), exact_options(8, {{"Sz", 0.5}})), std::invalid_argument);
}

TEST(ClassicDmrgValidation, WarmupConfigurationsMustFitTheChainAndTheSector) {
    auto options = exact_options(8, {{"Sz", 0.}});
    options.warmup_configuration = {0, 1, 0, 1, 0};
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument); // five states for six sites
    options.warmup_configuration = {0, 1, 0, 1, 0, 2};
    EXPECT_THROW(cd::solve(chain(), options), std::invalid_argument); // spin 1/2 has two states
    options.warmup_configuration = {0, 0, 0, 1, 0, 1};
    try {
        cd::solve(chain(), options);
        ADD_FAILURE() << "a configuration with Sz = 1 was accepted for the sector Sz = 0";
    } catch (std::invalid_argument const& error) {
        EXPECT_NE(std::string(error.what()).find("Sz = 1, but the sector asks for 0"), std::string::npos)
            << error.what();
    }
    cd::Evolution evolution;
    evolution.steps = 1;
    EXPECT_THROW(cd::evolve(chain(), chain(), evolution, options), std::invalid_argument);
    options.warmup_configuration = {0, 1, 0, 1, 0, 1};
    EXPECT_NO_THROW(cd::solve(chain(), options));
}

TEST(ClassicDmrgValidation, TermShapesTheEngineCannotRepresentAreRejected) {
    auto const site = mdl::spin(0.5);
    auto three = chain();
    three.add(mdl::TermSum::located(site, "Sz", 0) * mdl::TermSum::located(site, "Sz", 1)
              * mdl::TermSum::located(site, "Sz", 2) * std::complex<double>(0.1));
    EXPECT_THROW(cd::solve(three, exact_options(8, {{"Sz", 0.}})), std::invalid_argument);

    auto const fermions = mdl::spinless_fermion_model(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    auto const basis = mdl::physical_basis(fermions);
    auto const c = mdl::spinless_fermion();
    auto odd = mdl::compile_hamiltonian(fermions).op;
    odd.terms.push_back(mdl::LocalTerm{0.1, {mdl::LocalFactor{2, c->op("c").matrix}}, {0, 1}});
    EXPECT_THROW(cd::solve(odd, basis, exact_options(8, {{"N", 3.}})), std::invalid_argument);

    auto const hop = fermions.operator_decomposition(mdl::TermSum::located(c, "c_dag", 0)
                                                     * mdl::TermSum::located(c, "c", 1)).op;
    EXPECT_THROW(cd::solve(hop, basis, exact_options(8, {{"N", 3.}})), std::invalid_argument); // no h.c.

    auto wrong_string = hop;
    wrong_string.terms.front().jordan_wigner = {};
    EXPECT_THROW(cd::solve(wrong_string, basis, exact_options(8, {{"N", 3.}})), std::invalid_argument);
}

TEST(ClassicDmrgValidation, OperatorsMustConserveTheRequestedCharges) {
    auto const model = mdl::transverse_field_ising(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 0.7);
    EXPECT_THROW(cd::solve(model, exact_options(8, {{"Sz", 0.}})), std::invalid_argument); // not declared
    auto basis = mdl::physical_basis(model);
    basis.conserved_quantum_numbers = {"Sz"};
    EXPECT_THROW(cd::solve(mdl::compile_hamiltonian(model).op, basis, exact_options(8, {{"Sz", 0.}})),
                 std::invalid_argument);
}

TEST(ClassicDmrgValidation, MeasurementsMustLiveOnTheHamiltonianSites) {
    cd::Measurements requests;
    requests.terms["other"] = cd::site_terms(chain(8), "Sz");
    EXPECT_THROW(cd::solve(chain(6), exact_options(8, {{"Sz", 0.}}), requests), std::invalid_argument);
    requests.terms.clear();
    requests.terms[""] = cd::site_terms(chain(6), "Sz");
    EXPECT_THROW(cd::solve(chain(6), exact_options(8, {{"Sz", 0.}}), requests), std::invalid_argument);
}

TEST(ClassicDmrgValidation, RecursiveInvocationIsRefused) {
    auto const model = chain();
    EXPECT_THROW(cd::solve(model, exact_options(8, {{"Sz", 0.}}), {},
                           [&](cd::Iteration const&) { cd::solve(model, exact_options(8, {{"Sz", 0.}})); }),
                 std::logic_error);
    EXPECT_NO_THROW(cd::solve(model, exact_options(8, {{"Sz", 0.}})));
}
