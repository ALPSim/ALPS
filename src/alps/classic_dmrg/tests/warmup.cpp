// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Unevenly distributed ground states with Options::warmup_configuration on
// chains too long for ED, against exact references: fields alone (the
// ground state is the best configuration, found by sorting) and the XX
// chain in fields (free fermions, one tridiagonal matrix). The uniform
// warm-up leaves each of these several charges away with the default four
// sweeps (the 48-site wall at E = -8 instead of -24); the configuration's
// prefix charges place them before the sweeps start. A Hubbard potential
// step covers two charges at once and evolve() the preparation of its
// initial state.
#include "support.hpp"
#include <gtest/gtest.h>

#include <alps/models/compiled.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

using namespace classic_test;

namespace {

int const length = 48;

mdl::CompiledOperator sum(mdl::CompiledOperator a, mdl::CompiledOperator const& b) {
    a.terms.insert(a.terms.end(), b.terms.begin(), b.terms.end());
    return a;
}

// sum_i h_i S^z_i on the sites of `model`.
mdl::CompiledOperator fields(mdl::Hamiltonian const& model, std::vector<double> const& h) {
    auto op = cd::site_terms(model, "Sz");
    for (auto& term : op.terms) term.coefficient *= h.at(static_cast<std::size_t>(term.factors[0].site));
    return op;
}

// -1 left of `wall`, +1 from it on: the ground state is up there, down here.
std::vector<double> wall_fields(int wall) {
    std::vector<double> h(length);
    for (int i = 0; i < length; ++i) h[static_cast<std::size_t>(i)] = i < wall ? -1. : 1.;
    return h;
}

// Fields spread over (-1, 1) without a pattern (a Weyl sequence).
std::vector<double> irregular_fields() {
    std::vector<double> h(length);
    for (int i = 0; i < length; ++i) h[static_cast<std::size_t>(i)] = 2. * std::fmod(0.3247179572 * (i + 1), 1.) - 1.;
    return h;
}

// Spin up (local state 0) on the `up` sites with the lowest fields.
std::vector<std::size_t> best_configuration(std::vector<double> const& h, int up) {
    std::vector<std::size_t> order(h.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return h[a] < h[b]; });
    std::vector<std::size_t> configuration(h.size(), 1);
    for (int k = 0; k < up; ++k) configuration[order[static_cast<std::size_t>(k)]] = 0;
    return configuration;
}

// The XX chain J (S^x S^x + S^y S^y) in fields h, `up` spins up: free
// fermions hopping by J/2 in the potential h, less sum(h)/2.
double free_fermion_energy(std::vector<double> const& h, double J, int up) {
    int const n = static_cast<int>(h.size());
    Eigen::MatrixXd single = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) single(i, i) = h[static_cast<std::size_t>(i)];
    for (int i = 0; i + 1 < n; ++i) single(i, i + 1) = single(i + 1, i) = J / 2.;
    Eigen::VectorXd const levels = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(single).eigenvalues();
    return levels.head(up).sum() - std::accumulate(h.begin(), h.end(), 0.) / 2.;
}

cd::Options options_for(std::vector<std::size_t> configuration, int up) {
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(4, 32);
    options.sector = {{"Sz", up - length / 2.}};
    options.warmup_configuration = std::move(configuration);
    return options;
}

struct Case {
    double energy;
    std::vector<double> sz;
};

Case solve(std::vector<double> const& h, double J, int up) {
    auto const xx = mdl::xxz(lat::chain(length, 1.0, lat::Boundary::open), 1.0, 0.0);
    auto op = fields(xx, h);
    if (J != 0.) {
        auto hopping = mdl::compile_hamiltonian(xx).op;
        for (auto& term : hopping.terms) term.coefficient *= J;
        op = sum(op, hopping);
    }
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(xx, "Sz");
    auto const result = cd::solve(op, mdl::physical_basis(xx), options_for(best_configuration(h, up), up), requests);
    return {result.energies.at(0), result.measurements.at("Sz").values.at(0)};
}

} // namespace

TEST(ClassicDmrgWarmup, FieldOnlyWallsAreTheirConfiguration) {
    for (int wall : {length / 2, length / 4}) {
        auto const h = wall_fields(wall);
        auto const found = solve(h, 0., wall);
        EXPECT_NEAR(found.energy, -length / 2., 1e-10) << "wall at " << wall;
        ASSERT_EQ(found.sz.size(), static_cast<std::size_t>(length));
        for (int i = 0; i < length; ++i)
            EXPECT_NEAR(found.sz[static_cast<std::size_t>(i)], i < wall ? 0.5 : -0.5, 1e-10) << wall << " " << i;
    }
}

TEST(ClassicDmrgWarmup, IrregularFieldsReachTheBestConfiguration) {
    auto const h = irregular_fields();
    auto const configuration = best_configuration(h, length / 2);
    double exact = 0.;
    for (std::size_t i = 0; i < h.size(); ++i) exact += (configuration[i] == 0 ? 0.5 : -0.5) * h[i];
    EXPECT_NEAR(solve(h, 0., length / 2).energy, exact, 1e-10);
}

TEST(ClassicDmrgWarmup, WeakHoppingConvergesFromTheConfiguration) {
    // The configuration only places the charges; the sweeps then find the
    // quantum ground state (the uniform warm-up gives -18.56 for the wall).
    auto const wall = wall_fields(length / 4);
    EXPECT_NEAR(solve(wall, 0.1, length / 4).energy, free_fermion_energy(wall, 0.1, length / 4), 1e-8);
    auto const irregular = irregular_fields();
    EXPECT_NEAR(solve(irregular, 1e-3, length / 2).energy, free_fermion_energy(irregular, 1e-3, length / 2), 1e-8);
}

TEST(ClassicDmrgWarmup, HubbardPotentialStepSetsBothCharges) {
    // No hopping, U = 4 and -2 per particle on the left half: the eight
    // particles fill the left half singly (E = -16); the configuration puts
    // up spins (local state 1) on sites 0..3 and down spins (2) on 4..7.
    int const sites = 16;
    auto const model = mdl::hubbard(lat::chain(sites, 1.0, lat::Boundary::open), 0.0, 4.0);
    auto potential = cd::site_terms(model, "n");
    for (auto& term : potential.terms) term.coefficient *= term.factors[0].site < sites / 2 ? -2. : 0.;
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(4, 32);
    options.sector = {{"N", 8.}, {"Sz", 0.}};
    options.warmup_configuration.assign(static_cast<std::size_t>(sites), 0);
    for (int i = 0; i < sites / 2; ++i) options.warmup_configuration[static_cast<std::size_t>(i)] = i < sites / 4 ? 1 : 2;
    cd::Measurements requests;
    requests.terms["n"] = cd::site_terms(model, "n");
    auto const result = cd::solve(sum(mdl::compile_hamiltonian(model).op, potential), mdl::physical_basis(model), options,
                                  requests);
    EXPECT_NEAR(result.energies.at(0), -16., 1e-10);
    auto const& n = result.measurements.at("n").values.at(0);
    for (int i = 0; i < sites; ++i) EXPECT_NEAR(n[static_cast<std::size_t>(i)], i < sites / 2 ? 1. : 0., 1e-10) << i;
    EXPECT_EQ(result.options.warmup_configuration, options.warmup_configuration);
}

TEST(ClassicDmrgWarmup, EvolutionPreparesItsInitialStateFromTheConfiguration) {
    int const wall = length / 4;
    auto const xx = mdl::xxz(lat::chain(length, 1.0, lat::Boundary::open), 1.0, 0.0);
    cd::Evolution evolution;
    evolution.time_step = 0.05;
    evolution.steps = 1;
    evolution.states = 32;
    evolution.measure_steps = {0};
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(xx, "Sz");
    auto const result = cd::evolve(fields(xx, wall_fields(wall)), mdl::compile_hamiltonian(xx).op, mdl::physical_basis(xx),
                                   evolution, options_for(best_configuration(wall_fields(wall), wall), wall), requests);
    EXPECT_NEAR(result.initial.energies.at(0), -length / 2., 1e-10);
    auto const& sz = result.measurements.at("Sz").values.at(0);
    for (int i = 0; i < length; ++i) EXPECT_NEAR(sz[static_cast<std::size_t>(i)], i < wall ? 0.5 : -0.5, 1e-10) << i;
}
