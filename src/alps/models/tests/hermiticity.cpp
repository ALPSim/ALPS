// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/hermiticity.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
namespace mdl = alps::models;
using C = mdl::Complex;
constexpr auto strict = mdl::FermionConvention::strict;

mdl::TermSum op(mdl::SiteBasisPtr const& site, char const* name, int index = 0) {
    return mdl::TermSum::located(site, name, index);
}

mdl::SiteBasisPtr projectors() {
    mdl::SiteBasis site = *mdl::spin();
    site.operators["Pup"] = {"Pup", 2, {1., 0., 0., 0.}};
    site.operators["Pdown"] = {"Pdown", 2, {0., 0., 0., 1.}};
    site.operators["TinySp"] = {"TinySp", 2, {0., 0x1p-1000, 0., 0.}};
    return mdl::finalize_site(std::move(site));
}

mdl::TermSum long_product(mdl::SiteBasisPtr const& site, int count) {
    auto result = op(site, "Sz", 0);
    for (int i = 1; i < count; ++i) result = result * op(site, "Sz", i);
    return result;
}
} // namespace

TEST(Hermiticity, RejectsUnequalTransitionsAndImaginaryDiagonals) {
    auto site = mdl::spin();
    // Independent one-site matrices: [[0,1],[2,0]] and diag(i/2,-i/2).
    EXPECT_THROW(mdl::validate_hermiticity(op(site, "Sp") + 2. * op(site, "Sm"), {site}),
                 std::invalid_argument);
    EXPECT_THROW(mdl::validate_hermiticity(C(0., 1.) * op(site, "Sz"), {site}),
                 std::invalid_argument);
    auto const amplitude = C(.75, .25);
    EXPECT_NO_THROW(mdl::validate_hermiticity(
        amplitude * op(site, "Sp") + std::conj(amplitude) * op(site, "Sm"), {site}, strict, 0.));
}

TEST(Hermiticity, ExactCrossSupportCancellationIsRecognized) {
    auto site = projectors();
    // Pup + Pdown = I, so the first three terms cancel as a full operator.
    // The cancellation spans different tensor supports and is not termwise
    // Hermiticity. Its exact result is the explicitly Hermitian Sz(0)Sz(1).
    std::vector<mdl::TermSum> pieces{
        0x1p54 * op(site, "Sp", 0) * op(site, "Pup", 1),
        0x1p54 * op(site, "Sp", 0) * op(site, "Pdown", 1),
        -0x1p54 * op(site, "Sp", 0)};
    for (bool reverse : {false, true}) {
        if (reverse) std::reverse(pieces.begin(), pieces.end());
        auto sum = op(site, "Sz", 0) * op(site, "Sz", 1);
        for (auto const& piece : pieces) sum = sum + piece;
        EXPECT_NO_THROW(mdl::validate_hermiticity(sum, {site, site}, strict, 0.));
        EXPECT_THROW(mdl::validate_hermiticity(sum + op(site, "Sp", 0), {site, site}),
                     std::invalid_argument);
    }
}

TEST(Hermiticity, LargeCancellationNeitherLosesNorInventsAnAsymmetricRemainder) {
    auto site = mdl::spin();
    for (bool cancels : {false, true}) {
        std::vector<double> coefficients{0x1p54, 1., -0x1p54};
        if (cancels) coefficients.push_back(-1.);
        std::sort(coefficients.begin(), coefficients.end());
        do {
            auto sum = op(site, "Sz");
            for (double value : coefficients) sum = sum + value * op(site, "Sp");
            if (cancels) {
                EXPECT_NO_THROW(mdl::validate_hermiticity(sum, {site}, strict, 0.));
            } else {
                EXPECT_THROW(mdl::validate_hermiticity(sum, {site}), std::invalid_argument);
            }
        } while (std::next_permutation(coefficients.begin(), coefficients.end()));
    }
}

TEST(Hermiticity, LikeCoefficientExportCannotHideCrossSupportImaginaryScalar) {
    auto site = projectors();
    // Pup + Pdown = I. The full operator is iI; coalescing the first
    // two coefficients into a double before validation would erase it.
    auto terms = C(0., 0x1p54) * op(site, "I") + C(0., 1.) * op(site, "I")
        - C(0., 0x1p54) * op(site, "Pup") - C(0., 0x1p54) * op(site, "Pdown");
    EXPECT_THROW(mdl::validate_hermiticity(terms, {site}), std::invalid_argument);
    EXPECT_NO_THROW(mdl::validate_hermiticity(terms - C(0., 1.) * op(site, "I"),
                                            {site}, strict, 0.));
}

TEST(Hermiticity, NoAbsoluteAmplitudeFloorOrCoefficientPlacementDependence) {
    auto site = projectors();
    // Both expressions are exactly 2^-900 Sp. Squared double magnitudes
    // underflow, but a finite unpaired transition is still non-Hermitian.
    auto const direct = 0x1p-900 * op(site, "Sp");
    auto const factored = 0x1p100 * op(site, "TinySp");
    for (auto const& transition : {direct, factored}) {
        EXPECT_THROW(mdl::validate_hermiticity(transition, {site}), std::invalid_argument);
        EXPECT_NO_THROW(mdl::validate_hermiticity(transition + transition.dag(), {site}, strict, 0.));
    }
}

TEST(Hermiticity, RealIdentityShiftsCannotMaskAsymmetricTermsOrImaginaryScalars) {
    auto site = mdl::spin();
    for (double shift : {0., 0x1p50, 0x1p900}) {
        auto const identity = shift * op(site, "I");
        EXPECT_NO_THROW(mdl::validate_hermiticity(identity, {site}, strict, 0.));
        EXPECT_THROW(mdl::validate_hermiticity(identity + 0x1p-500 * op(site, "Sp"), {site}),
                     std::invalid_argument);
        // With no nonidentity scale, even an arbitrarily small imaginary
        // scalar is rejected; the real energy origin is irrelevant.
        EXPECT_THROW(mdl::validate_hermiticity(identity + C(0., 0x1p-500) * op(site, "I"), {site}),
                     std::invalid_argument);
        EXPECT_THROW(mdl::validate_hermiticity(identity + op(site, "Sz")
            + C(0., 1e-9) * op(site, "I"), {site}), std::invalid_argument);
        EXPECT_NO_THROW(mdl::validate_hermiticity(identity + op(site, "Sz")
            + C(0., 1e-12) * op(site, "I"), {site}));
    }
}

TEST(Hermiticity, DeclaredPairRelativeToleranceHasAnInclusiveBoundary) {
    auto site = mdl::spin();
    double const tolerance = 0x1p-20;
    double const inside = 1. - tolerance;
    double const outside = std::nextafter(inside, 0.);
    for (double scale : {0x1p-900, 1., 0x1p900}) {
        // The two offdiagonal magnitudes are scale and scale*(1-tol).
        // The residual equals tol*scale exactly, with no unrelated scale.
        auto const boundary = scale * op(site, "Sp") + (scale * inside) * op(site, "Sm");
        auto const beyond = scale * op(site, "Sp") + (scale * outside) * op(site, "Sm");
        EXPECT_NO_THROW(mdl::validate_hermiticity(boundary, {site}, strict, tolerance));
        EXPECT_THROW(mdl::validate_hermiticity(beyond, {site}, strict, tolerance), std::invalid_argument);
        EXPECT_THROW(mdl::validate_hermiticity(boundary, {site}, strict, 0.), std::invalid_argument);
    }
    EXPECT_NO_THROW(mdl::validate_hermiticity(op(site, "Sp") + (1. - 0x1p-38) * op(site, "Sm"), {site}));
    EXPECT_THROW(mdl::validate_hermiticity(op(site, "Sp") + (1. - 0x1p-35) * op(site, "Sm"), {site}),
                 std::invalid_argument);
}

TEST(Hermiticity, UsesComplexPairMagnitudeRatherThanIndependentComponentCutoffs) {
    auto site = mdl::spin();
    C const amplitude(3., 4.);
    double const tolerance = 0x1p-20;
    auto const near = amplitude * op(site, "Sp")
        + std::conj(amplitude * (1. - tolerance)) * op(site, "Sm");
    auto const far = amplitude * op(site, "Sp")
        + std::conj(amplitude * (1. - 2. * tolerance)) * op(site, "Sm");
    EXPECT_NO_THROW(mdl::validate_hermiticity(near, {site}, strict, tolerance));
    EXPECT_THROW(mdl::validate_hermiticity(far, {site}, strict, tolerance), std::invalid_argument);
}

TEST(Hermiticity, RespectsFermionConventionAndReversedFactorOrder) {
    auto site = mdl::spinless_fermion();
    C const amplitude(.75, .25);
    auto const forward = amplitude * op(site, "c_dag", 2) * op(site, "c", 0);
    auto const adjoint = std::conj(amplitude) * op(site, "c_dag", 0) * op(site, "c", 2);
    auto const reversed = amplitude * op(site, "c", 0) * op(site, "c_dag", 2);
    for (auto convention : {strict, mdl::FermionConvention::lattice}) {
        SCOPED_TRACE(convention == strict ? "strict" : "lattice");
        EXPECT_NO_THROW(mdl::validate_hermiticity(forward + adjoint, {site, site, site}, convention, 0.));
        EXPECT_THROW(mdl::validate_hermiticity(forward + 2. * adjoint, {site, site, site}, convention),
                     std::invalid_argument);
        double const sign = convention == strict ? 1. : -1.;
        // Distinct-site CAR gives c0 cdag2 = -cdag2 c0 in strict
        // convention. The lattice tensor convention has no swap sign.
        EXPECT_NO_THROW(mdl::validate_hermiticity(forward + sign * reversed, {site, site, site}, convention, 0.));
        EXPECT_THROW(mdl::validate_hermiticity(forward - sign * reversed, {site, site, site}, convention),
                     std::invalid_argument);
    }
}

TEST(Hermiticity, DoesNotImposeConservationOrRestrictGenericCompilation) {
    auto site = mdl::spin();
    mdl::Hamiltonian declared(alps::lattice::chain(2), site, strict, {"Sz"});
    declared.add(op(site, "Sx", 0));
    EXPECT_NO_THROW(mdl::validate_hermiticity(declared));
    EXPECT_THROW(mdl::compile_hamiltonian(declared), std::invalid_argument);

    auto const arbitrary = op(site, "Sp") + 2. * op(site, "Sm");
    auto const before = mdl::compile_hamiltonian(arbitrary, {site});
    EXPECT_THROW(mdl::validate_hermiticity(arbitrary, {site}), std::invalid_argument);
    auto const after = mdl::compile_hamiltonian(arbitrary, {site});
    ASSERT_EQ(after.terms.size(), 2);
    ASSERT_EQ(before.terms.size(), 2);
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_EQ(after.terms[i].coefficient, before.terms[i].coefficient);
        EXPECT_EQ(after.terms[i].jordan_wigner, before.terms[i].jordan_wigner);
        ASSERT_EQ(after.terms[i].factors.size(), 1);
        ASSERT_EQ(before.terms[i].factors.size(), 1);
        EXPECT_EQ(after.terms[i].factors[0].site, before.terms[i].factors[0].site);
        EXPECT_EQ(after.terms[i].factors[0].matrix, before.terms[i].factors[0].matrix);
    }
    EXPECT_EQ(after.terms[0].coefficient, C(1.));
    EXPECT_EQ(after.terms[1].coefficient, C(2.));
    EXPECT_EQ(after.terms[0].factors[0].matrix, (mdl::Matrix{0., 1., 0., 0.}));
    EXPECT_EQ(after.terms[1].factors[0].matrix, (mdl::Matrix{0., 0., 1., 0.}));
}

TEST(Hermiticity, LongHermitianProductsDoNotExpandTheGlobalBasis) {
    auto site = mdl::spin();
    auto const product = long_product(site, 64);
    // Sixty-four Sz factors have 2^64 independent-word contributions,
    // but their Hermiticity is structural. A tiny fallback budget ensures
    // this test cannot pass by enumerating them.
    EXPECT_NO_THROW(mdl::validate_hermiticity(product, std::vector<mdl::SiteBasisPtr>(64, site), strict, 0., 1));
    auto fermion = mdl::spinless_fermion();
    auto const hopping = C(.75, .25) * op(fermion, "c_dag", 63) * op(fermion, "c", 0);
    for (auto convention : {strict, mdl::FermionConvention::lattice}) {
        EXPECT_NO_THROW(mdl::validate_hermiticity(hopping.add_hc(),
            std::vector<mdl::SiteBasisPtr>(64, fermion), convention, 0., 1));
        // The asymmetric counterpart needs the general coefficient check.
        // Its long JW string must still not consume an exponential budget.
        EXPECT_THROW(mdl::validate_hermiticity(hopping + 2. * hopping.dag(),
            std::vector<mdl::SiteBasisPtr>(64, fermion), convention, 0., 16), std::invalid_argument);
    }
}

TEST(Hermiticity, AmbiguousExpansionIsBoundedAndDoesNotChangeTheOperator) {
    auto site = projectors();
    auto const hermitian = op(site, "Sp", 0) * op(site, "Pup", 1)
        + op(site, "Sp", 0) * op(site, "Pdown", 1) + op(site, "Sm", 0);
    // Algebraically Sp+Sm. The explicit small resource limit must not
    // silently replace its cross-support cancellation with a projection.
    EXPECT_THROW(mdl::validate_hermiticity(hermitian, {site, site}, strict, 0., 1), std::length_error);
    EXPECT_NO_THROW(mdl::validate_hermiticity(hermitian, {site, site}, strict, 0., 128));
    auto const long_ambiguous = long_product(site, 64) + op(site, "Sp", 0);
    // A future structural rejection is equally correct. This must not be
    // accepted, allocate the global tensor expansion, or fail with bad_alloc.
    try {
        mdl::validate_hermiticity(long_ambiguous,
            std::vector<mdl::SiteBasisPtr>(64, site), strict, 0., 32);
        FAIL() << "asymmetric long operator was accepted";
    } catch (std::length_error const&) {
    } catch (std::invalid_argument const&) {
    }
}

TEST(Hermiticity, LocalSixtyFourSiteModelsRetainBothOverloads) {
    auto lattice = alps::lattice::chain(64, 1., alps::lattice::Boundary::open);
    auto const heisenberg = mdl::heisenberg(lattice);
    EXPECT_NO_THROW(mdl::validate_hermiticity(heisenberg, 0.));
    auto site = mdl::spin();
    mdl::Hamiltonian xy(lattice, site, strict, {"Sz"});
    for (int i = 0; i < 63; ++i)
        xy.add(op(site, "Sx", i) * op(site, "Sx", i + 1)
             + op(site, "Sy", i) * op(site, "Sy", i + 1));
    EXPECT_NO_THROW(mdl::validate_hermiticity(xy, 0.));
    EXPECT_NO_THROW(mdl::validate_hermiticity(xy.expand(),
        std::vector<mdl::SiteBasisPtr>(64, site), strict, 0.));
}

TEST(Hermiticity, RejectsInvalidSettingsAndMalformedNumericalInputs) {
    auto site = mdl::spin();
    auto const valid = op(site, "Sz");
    for (double tolerance : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        EXPECT_THROW(mdl::validate_hermiticity(valid, {site}, strict, tolerance), std::invalid_argument);
    EXPECT_THROW(mdl::validate_hermiticity(valid, {site}, strict, 0., 0), std::invalid_argument);
    EXPECT_THROW(mdl::validate_hermiticity(valid, {site}, static_cast<mdl::FermionConvention>(99)),
                 std::invalid_argument);
    EXPECT_THROW(mdl::validate_hermiticity(valid, {nullptr}), std::invalid_argument);
    EXPECT_THROW(mdl::validate_hermiticity(op(site, "Sz", 1), {site}), std::invalid_argument);
    EXPECT_THROW(mdl::validate_hermiticity(std::numeric_limits<double>::infinity() * valid, {site}),
                 std::invalid_argument);
}
