// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The matrix product states of Options::keep_mps against exact
// diagonalization: contracted to a vector on the product basis, every
// targeted state must be the ED eigenvector up to a phase, with the
// canonical form the header states. Spin-1/2 and spin-1 chains, the Hubbard
// chain, frustrated spinless fermions (which lock the fermion ordering of the
// amplitudes), complex hopping phases, excited targets with and without
// stored products (the two closing configurations) and odd lengths.
#include "support.hpp"
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <variant>
#include <vector>

using namespace classic_test;

namespace {

using Complex = std::complex<double>;

std::vector<Complex> entries(cd::MpsTensor const& site) {
    if (!site.complex.empty()) return site.complex;
    return std::vector<Complex>(site.real.begin(), site.real.end());
}

// The state on the whole product basis, site 0 most significant.
std::vector<Complex> contract(cd::MatrixProductState const& mps) {
    std::vector<Complex> carry{1.}; // [configuration][bond]
    std::size_t bond = 1;
    for (auto const& site : mps.sites) {
        EXPECT_EQ(site.left, bond);
        auto const t = entries(site);
        EXPECT_EQ(t.size(), site.left * site.local * site.right);
        std::size_t const configurations = carry.size() / bond;
        std::vector<Complex> next(configurations * site.local * site.right);
        for (std::size_t c = 0; c < configurations; ++c)
            for (std::size_t l = 0; l < site.left; ++l) {
                Complex const x = carry[c * bond + l];
                if (x == Complex{}) continue;
                for (std::size_t s = 0; s < site.local; ++s)
                    for (std::size_t r = 0; r < site.right; ++r)
                        next[(c * site.local + s) * site.right + r] += x * t[(l * site.local + s) * site.right + r];
            }
        carry = std::move(next);
        bond = site.right;
    }
    EXPECT_EQ(bond, 1u);
    return carry;
}

// sum_{l, s} conj(A[l, s, x]) A[l, s, y] = delta_xy (left) or the mirror (right).
double orthonormality_error(cd::MpsTensor const& site, bool left) {
    auto const t = entries(site);
    std::size_t const k = left ? site.right : site.left;
    double worst = 0.;
    for (std::size_t x = 0; x < k; ++x)
        for (std::size_t y = 0; y < k; ++y) {
            Complex sum{};
            for (std::size_t a = 0; a < (left ? site.left : site.right); ++a)
                for (std::size_t s = 0; s < site.local; ++s) {
                    auto at = [&](std::size_t bond) {
                        return left ? t[(a * site.local + s) * site.right + bond]
                                    : t[(bond * site.local + s) * site.right + a];
                    };
                    sum += std::conj(at(x)) * at(y);
                }
            worst = std::max(worst, std::abs(sum - (x == y ? 1. : 0.)));
        }
    return worst;
}

std::vector<Complex> ed_vector(alps::ed::ProductEigenResult const& states, std::size_t column) {
    std::size_t const dim = states.problem.dimension();
    return std::visit([&](auto const& values) {
        std::vector<Complex> out;
        for (std::size_t i = 0; i < dim; ++i) out.emplace_back(values[column * dim + i]);
        return out;
    }, states.eigenvectors);
}

// Every targeted state against ED's eigenvector of the same rank.
void expect_exact_states(mdl::Hamiltonian const& model, cd::Options options, Sector const& sector,
                         std::size_t center, cd::Measurements const& requests = {}) {
    options.keep_mps = true;
    auto const result = cd::solve(model, options, requests);
    auto const states = ed_states(model, sector, static_cast<int>(options.num_eigenvalues));
    ASSERT_EQ(result.mps.size(), options.num_eigenvalues);
    auto const& dims = states.problem.site_dims();
    auto const& configurations = states.problem.states();
    std::size_t const length = dims.size();
    for (std::size_t target = 0; target < result.mps.size(); ++target) {
        SCOPED_TRACE(target);
        auto const& mps = result.mps[target];
        ASSERT_EQ(mps.sites.size(), length);
        EXPECT_EQ(mps.center, center);
        for (std::size_t i = 0; i < length; ++i) {
            EXPECT_EQ(mps.sites[i].local, static_cast<std::size_t>(dims[i])) << i;
            EXPECT_EQ(mps.sites[i].complex.empty(), options.arithmetic == cd::Arithmetic::real) << i;
            if (i != mps.center) EXPECT_LT(orthonormality_error(mps.sites[i], i < mps.center), 1e-10) << i;
        }
        auto const dense = contract(mps);
        auto const reference = ed_vector(states, target);
        // The overlap with the ED vector, on its sector states, and the weight outside.
        Complex overlap{};
        double norm = 0., inside = 0.;
        for (auto const& x : dense) norm += std::norm(x);
        for (std::size_t k = 0; k < reference.size(); ++k) {
            std::size_t index = 0;
            for (std::size_t i = 0; i < length; ++i)
                index = index * static_cast<std::size_t>(dims[i]) + static_cast<std::size_t>(configurations[k * length + i]);
            overlap += std::conj(reference[k]) * dense[index];
            inside += std::norm(dense[index]);
        }
        EXPECT_NEAR(norm, 1., 1e-10);
        EXPECT_NEAR(inside, 1., 1e-10);
        EXPECT_NEAR(std::abs(overlap), 1., 1e-8);
        EXPECT_NEAR(result.energies[target], states.eigenvalues[target], 1e-8);
    }
}

} // namespace

TEST(ClassicDmrgMps, HeisenbergChainIsTheExactGroundState) {
    auto const model = mdl::heisenberg(lat::chain(10, 1.0, lat::Boundary::open), 1.0);
    expect_exact_states(model, exact_options(64, {{"Sz", 0.}}), {{"Sz", 0.}}, 5);
}

TEST(ClassicDmrgMps, OddLengthsCloseAtTheMiddle) {
    auto const model = mdl::heisenberg(lat::chain(9, 1.0, lat::Boundary::open), 1.0, 0.1);
    expect_exact_states(model, exact_options(64, {{"Sz", 0.5}}), {{"Sz", 0.5}}, 4);
}

TEST(ClassicDmrgMps, FourSitesAreTheSmallestChain) {
    auto const model = mdl::heisenberg(lat::chain(4, 1.0, lat::Boundary::open), 1.0);
    expect_exact_states(model, exact_options(16, {{"Sz", 0.}}), {{"Sz", 0.}}, 2);
}

TEST(ClassicDmrgMps, SpinOneSitesMapToTheirBasisOrder) {
    // The engine sorts each site basis by its charges; the tensors use the
    // PhysicalBasis order of the model's sites.
    auto const model = mdl::heisenberg_d(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 0.3, 0.2, 1.0);
    expect_exact_states(model, exact_options(81, {{"Sz", 1.}}), {{"Sz", 1.}}, 3);
}

TEST(ClassicDmrgMps, UnrestrictedStatesHaveNoSectorSubspaces) {
    auto const model = mdl::transverse_field_ising(lat::chain(8, 1.0, lat::Boundary::open), 1.0, 0.7);
    expect_exact_states(model, exact_options(32), {}, 4);
}

TEST(ClassicDmrgMps, HubbardAmplitudesFollowTheAlpsFermionOrder) {
    auto const model = mdl::hubbard(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 4.0, 0.2);
    Sector const sector{{"N", 5.}, {"Sz", 0.5}};
    expect_exact_states(model, exact_options(128, sector), sector, 3);
}

TEST(ClassicDmrgMps, FrustratedSpinlessFermionsKeepTheirSigns) {
    // Next-nearest-neighbour hopping closes triangles, so amplitudes in a
    // different fermion order would not be an eigenvector of this model.
    auto const model = mdl::spinless_fermion_model(lat::nnn_chain(8, 1.0, lat::Boundary::open), 1.0, 0.3, 0.4);
    expect_exact_states(model, exact_options(64, {{"N", 3.}}), {{"N", 3.}}, 4);
}

TEST(ClassicDmrgMps, ComplexHoppingPhasesGiveComplexTensors) {
    auto const site = mdl::spinless_fermion();
    mdl::Hamiltonian model(lat::nnn_chain(8, 1.0, lat::Boundary::open), site, mdl::FermionConvention::strict, {"N"});
    auto hop = [site](int i, int j) {
        return mdl::TermSum::located(site, "c_dag", i) * mdl::TermSum::located(site, "c", j);
    };
    model.add_bond(-std::polar(1.0, 0.4), hop, 0, true);
    model.add_bond(-0.5, hop, 1, true);
    auto options = exact_options(64, {{"N", 3.}});
    options.arithmetic = cd::Arithmetic::complex;
    expect_exact_states(model, options, {{"N", 3.}}, 4);
}

TEST(ClassicDmrgMps, ExcitedTargetsShareTheBlocks) {
    auto const model = mdl::heisenberg(lat::chain(8, 1.0, lat::Boundary::open), 1.0);
    expect_exact_states(model, exact_options(64, {{"Sz", 0.}}, 3), {{"Sz", 0.}}, 4);
}

TEST(ClassicDmrgMps, StoredProductsCloseAtTheRightEnd) {
    // Several targets with a two-site measurement store products and sweep
    // the whole chain, so the closing configuration is at the right end.
    auto const model = mdl::heisenberg(lat::chain(8, 1.0, lat::Boundary::open), 1.0);
    cd::Measurements requests;
    requests.terms["SzSz"] = cd::pair_terms(model, "Sz", "Sz");
    expect_exact_states(model, exact_options(64, {{"Sz", 0.}}, 2), {{"Sz", 0.}}, 6, requests);
}

TEST(ClassicDmrgMps, StoredProductsWithFermions) {
    auto const model = mdl::hubbard(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 2.0);
    Sector const sector{{"N", 6.}, {"Sz", 0.}};
    cd::Measurements requests;
    requests.terms["nn"] = cd::pair_terms(model, "n", "n");
    expect_exact_states(model, exact_options(128, sector, 2), sector, 4, requests);
}

TEST(ClassicDmrgMps, TruncatedStatesKeepTheirCanonicalForm) {
    auto const model = mdl::heisenberg(lat::chain(16, 1.0, lat::Boundary::open), 1.0);
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(3, 12);
    options.warmup_states = 12;
    options.sector = {{"Sz", 0.}};
    options.keep_mps = true;
    auto const result = cd::solve(model, options);
    ASSERT_EQ(result.mps.size(), 1u);
    auto const& mps = result.mps.front();
    EXPECT_EQ(mps.center, 8u);
    for (std::size_t i = 0; i < mps.sites.size(); ++i) {
        EXPECT_LE(std::max(mps.sites[i].left, mps.sites[i].right), 2u * 12u) << i;
        if (i != mps.center) EXPECT_LT(orthonormality_error(mps.sites[i], i < mps.center), 1e-10) << i;
    }
    double norm = 0.;
    for (auto const& x : contract(mps)) norm += std::norm(x);
    EXPECT_NEAR(norm, 1., 1e-10);
}

TEST(ClassicDmrgMps, OffByDefault) {
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    EXPECT_TRUE(cd::solve(model, exact_options(16, {{"Sz", 0.}})).mps.empty());
}
