// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ground and excited state energies of the classic DMRG solver against exact
// diagonalization: spin chains and ladders, a spin-1 chain, hard-core bosons,
// the Hubbard chain, frustrated spinless fermions (a non-bipartite hopping
// graph locks the Jordan-Wigner conventions) and complex hopping phases.
#include "support.hpp"
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <stdexcept>

using namespace classic_test;

TEST(ClassicDmrgEnergies, HeisenbergOpenChainMatchesExactDiagonalization) {
    auto const model = mdl::heisenberg(lat::chain(8, 1.0, lat::Boundary::open), 1.0);
    auto const result = cd::solve(model, exact_options(64, {{"Sz", 0.}}));
    ASSERT_EQ(result.energies.size(), 1u);
    EXPECT_NEAR(result.energies[0], ed_energies(model, {{"Sz", 0.}})[0], 1e-9);
    EXPECT_LT(result.final_truncation_error, 1e-12);
    EXPECT_EQ(result.options.sweeps, cd::fixed_sweeps(4, 64));
    EXPECT_EQ(result.options.sector.at("Sz"), 0.);
    EXPECT_GT(result.options.lanczos_tolerance, 0.);
}

TEST(ClassicDmrgEnergies, HeisenbergOpenLadderMatchesExactDiagonalization) {
    auto const model = mdl::heisenberg(lat::ladder(4, 1.0, lat::Boundary::open), 1.0);
    auto const result = cd::solve(model, exact_options(64, {{"Sz", 0.}}));
    EXPECT_NEAR(result.energies[0], ed_energies(model, {{"Sz", 0.}})[0], 1e-9);
}

TEST(ClassicDmrgEnergies, SpinOneChainWithSingleIonAnisotropy) {
    auto const model = mdl::heisenberg_d(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 0.3, 0.0, 1.0);
    auto const result = cd::solve(model, exact_options(81, {{"Sz", 0.}}));
    EXPECT_NEAR(result.energies[0], ed_energies(model, {{"Sz", 0.}})[0], 1e-9);
}

TEST(ClassicDmrgEnergies, SectorTargetsSelectTheTotalCharge) {
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    for (double sz : {0., 1., 2.}) {
        SCOPED_TRACE(sz);
        auto const result = cd::solve(model, exact_options(32, {{"Sz", sz}}));
        EXPECT_NEAR(result.energies[0], ed_energies(model, {{"Sz", sz}})[0], 1e-9);
    }
}

TEST(ClassicDmrgEnergies, FieldBreaksTheSectorSymmetry) {
    // E(Sz = +1) != E(Sz = -1) in a field: the targeted total has a sign.
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 0.3);
    for (double sz : {1., -1.}) {
        SCOPED_TRACE(sz);
        auto const result = cd::solve(model, exact_options(32, {{"Sz", sz}}));
        EXPECT_NEAR(result.energies[0], ed_energies(model, {{"Sz", sz}})[0], 1e-9);
    }
}

TEST(ClassicDmrgEnergies, UnrestrictedSolveFindsTheGlobalGroundState) {
    auto const model = mdl::transverse_field_ising(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 0.7);
    auto const result = cd::solve(model, exact_options(32));
    EXPECT_NEAR(result.energies[0], ed_energies(model, {})[0], 1e-9);
}

TEST(ClassicDmrgEnergies, ExcitedStatesAreTargetedTogether) {
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    auto const result = cd::solve(model, exact_options(64, {{"Sz", 0.}}, 3));
    auto const reference = ed_energies(model, {{"Sz", 0.}}, 3);
    ASSERT_EQ(result.energies.size(), 3u);
    for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(result.energies[i], reference[i], 1e-8) << i;
}

TEST(ClassicDmrgEnergies, HardcoreBosonsConserveTheParticleNumber) {
    auto const model = mdl::hardcore_boson_model(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 0.0, 0.5);
    auto const result = cd::solve(model, exact_options(32, {{"N", 3.}}));
    EXPECT_NEAR(result.energies[0], ed_energies(model, {{"N", 3.}})[0], 1e-9);
}

TEST(ClassicDmrgEnergies, HubbardOpenChainMatchesExactDiagonalization) {
    auto const model = mdl::hubbard(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 4.0);
    Sector const sector{{"N", 6.}, {"Sz", 0.}};
    auto const result = cd::solve(model, exact_options(128, sector));
    EXPECT_NEAR(result.energies[0], ed_energies(model, sector)[0], 1e-8);
}

TEST(ClassicDmrgEnergies, HubbardAwayFromHalfFillingWithPolarization) {
    auto const model = mdl::hubbard(lat::chain(5, 1.0, lat::Boundary::open), 1.0, 2.0, 0.3);
    Sector const sector{{"N", 3.}, {"Sz", 0.5}};
    auto const result = cd::solve(model, exact_options(128, sector));
    EXPECT_NEAR(result.energies[0], ed_energies(model, sector)[0], 1e-8);
}

TEST(ClassicDmrgEnergies, FrustratedSpinlessFermionsLockTheExchangeSigns) {
    // Next-nearest-neighbour hopping closes triangles, so a wrong fermion
    // sign convention changes the spectrum.
    auto const model = mdl::spinless_fermion_model(lat::nnn_chain(6, 1.0, lat::Boundary::open), 1.0, 0.0, 0.4);
    auto const result = cd::solve(model, exact_options(64, {{"N", 3.}}));
    EXPECT_NEAR(result.energies[0], ed_energies(model, {{"N", 3.}})[0], 1e-9);
}

TEST(ClassicDmrgEnergies, ComplexHoppingPhasesNeedComplexArithmetic) {
    using namespace std::complex_literals;
    auto const site = mdl::spinless_fermion();
    mdl::Hamiltonian model(lat::nnn_chain(6, 1.0, lat::Boundary::open), site, mdl::FermionConvention::strict, {"N"});
    std::complex<double> const t = std::polar(1.0, 0.4);
    model.add_bond(-t, [site](int i, int j) { return mdl::TermSum::located(site, "c_dag", i) * mdl::TermSum::located(site, "c", j); }, 0, true);
    model.add_bond(-0.5, [site](int i, int j) { return mdl::TermSum::located(site, "c_dag", i) * mdl::TermSum::located(site, "c", j); }, 1, true);
    auto options = exact_options(64, {{"N", 3.}});
    EXPECT_THROW(cd::solve(model, options), std::invalid_argument);
    options.arithmetic = cd::Arithmetic::complex;
    auto const result = cd::solve(model, options);
    EXPECT_NEAR(result.energies[0], ed_energies(model, {{"N", 3.}})[0], 1e-9);
    EXPECT_EQ(result.options.arithmetic, cd::Arithmetic::complex);
}

TEST(ClassicDmrgEnergies, TruncatedSweepsConvergeFromAbove) {
    auto const model = mdl::heisenberg(lat::chain(12, 1.0, lat::Boundary::open), 1.0);
    cd::Options options;
    options.sweeps = cd::ramped_sweeps(3, 24);
    options.warmup_states = 8;
    options.sector = {{"Sz", 0.}};
    auto const result = cd::solve(model, options);
    double const exact = ed_energies(model, {{"Sz", 0.}})[0];
    EXPECT_GT(result.energies[0], exact - 1e-9);
    EXPECT_NEAR(result.energies[0], exact, 1e-3);
    EXPECT_GT(result.final_truncation_error, 0.);
    EXPECT_EQ(result.options.sweeps.size(), 3u);
}

TEST(ClassicDmrgEnergies, IterationsAreObservedAndCanAbortTheSolve) {
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    std::size_t observed = 0;
    auto const result = cd::solve(model, exact_options(16, {{"Sz", 0.}}), {}, [&](cd::Iteration const& row) {
        ++observed;
        EXPECT_TRUE(std::isfinite(row.energy));
        EXPECT_GE(row.truncation_error, -1e-12); // 1 - sum of kept weights, to rounding
    });
    EXPECT_EQ(observed, result.iterations.size());
    EXPECT_GT(observed, 0u);
    bool both_directions = false;
    for (auto const& row : result.iterations)
        both_directions = both_directions || row.direction != result.iterations.front().direction;
    EXPECT_TRUE(both_directions);
    struct Abort {};
    EXPECT_THROW(cd::solve(model, exact_options(16, {{"Sz", 0.}}), {}, [](cd::Iteration const&) { throw Abort{}; }), Abort);
    // The engine state is restored after the aborted solve.
    auto const again = cd::solve(model, exact_options(16, {{"Sz", 0.}}));
    EXPECT_NEAR(again.energies[0], result.energies[0], 1e-10);
}
