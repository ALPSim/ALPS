// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/term_accumulation.cpp.
// Adaptation: <alps/ed/...> and <alps/models/...> include names; the operator
// action is read from the canonical CSR on the full product basis instead of
// alpsx's projected-basis application.
#include <alps/ed/product_basis.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/models/hamiltonian.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <complex>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace mdl = alps::models;
using C = mdl::Complex;
using Input = alps::models::CompiledModel;
constexpr auto strict = mdl::FermionConvention::strict;
constexpr double A = 0x1p54;
constexpr double B = 0x1p1023;

enum class Route { hamiltonian, local };

mdl::TermSum op(mdl::SiteBasisPtr const& site, char const* name, int index = 0) {
    return mdl::TermSum::located(site, name, index);
}

Input lower(Route route, mdl::TermSum const& terms, mdl::SiteBasisPtr const& site,
            int count = 1, std::vector<std::string> charges = {},
            mdl::FermionConvention convention = strict) {
    mdl::Hamiltonian model(alps::lattice::chain(count, 1., alps::lattice::Boundary::open),
                           site, convention, std::move(charges));
    model.add(terms);
    return route == Route::hamiltonian ? mdl::compile_hamiltonian(model)
                                        : model.to_local_decomposition();
}

// Exercise an installed public consumer, rather than imposing a descriptor
// term count/order or a test-local coefficient summation policy. With no
// charge targets, this basis is the complete product basis, with site zero
// most significant; the canonical CSR supplies every column. Expected
// entries in the tests are elementary matrices/Fock-state amplitudes, never
// another compiler's output or a multiplication of the input factors.
void expect_action(Input const& input, mdl::SiteBasisPtr const& site,
                   std::vector<C> const& expected) {
    auto const basis = alps::ed::enumerate_product_basis(input.op.site_dims);
    auto const size = basis.dimension;
    ASSERT_EQ(expected.size(), size * size);
    std::vector<std::vector<double>> const parity(input.op.site_dims.size(), site->parity_operator);
    std::vector<std::int64_t> const states(basis.states.begin(), basis.states.end());
    auto const matrix = alps::ed::build_csr_from_local<C>(input.op, states, size, parity);
    std::vector<C> actual(size * size);
    for (std::size_t row = 0; row < size; ++row)
        for (auto p = matrix.offsets[row]; p < matrix.offsets[row + 1]; ++p)
            actual[row * size + static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(p)])]
                = matrix.values[static_cast<std::size_t>(p)];
    for (std::size_t row = 0; row < size; ++row)
        for (std::size_t column = 0; column < size; ++column) {
            SCOPED_TRACE(row);
            SCOPED_TRACE(column);
            // All oracle values are exactly representable binary numbers.
            EXPECT_EQ(actual[row * size + column], expected[row * size + column]);
        }
}

mdl::SiteBasisPtr product_site(C remainder = 1.) {
    mdl::SiteBasis site;
    site.dim = 3;
    site.quantum_numbers["Q"] = {0., 1., 2.};
    site.basis_labels = {"zero", "one", "two"};
    // X Y = remainder E01, since row 0 dot column 1 is A+r-A.
    site.operators["X"] = {"X", 3, {C(A, A), remainder, C(-A, -A),
                                       0., 0., 0., 0., 0., 0.}};
    site.operators["Y"] = {"Y", 3, {0., 1., 0., 0., 1., 0., 0., 1., 0.}};
    site.operators["Z"] = {"Z", 3, {0., 0., 0., 1., 0., 0., 0., 0., 0.}};
    return mdl::finalize_site(std::move(site));
}

void expect_unchanged(mdl::TermSum const& actual, mdl::TermSum const& before) {
    ASSERT_EQ(actual.terms().size(), before.terms().size());
    for (std::size_t t = 0; t < actual.terms().size(); ++t) {
        auto const& a = actual.terms()[t];
        auto const& b = before.terms()[t];
        EXPECT_EQ(a.coeff, b.coeff);
        EXPECT_EQ(a.sign, b.sign);
        ASSERT_EQ(a.factors.size(), b.factors.size());
        for (std::size_t f = 0; f < a.factors.size(); ++f) {
            EXPECT_EQ(a.factors[f].site, b.factors[f].site);
            EXPECT_EQ(a.factors[f].basis, b.factors[f].basis);
            EXPECT_EQ(a.factors[f].op.name, b.factors[f].op.name);
            EXPECT_EQ(a.factors[f].op.dim, b.factors[f].op.dim);
            EXPECT_EQ(a.factors[f].op.matrix, b.factors[f].op.matrix);
            EXPECT_EQ(a.factors[f].op.parity, b.factors[f].op.parity);
            EXPECT_EQ(a.factors[f].op.qn_deltas, b.factors[f].op.qn_deltas);
        }
    }
}

class TermAccumulation : public ::testing::TestWithParam<Route> {};
} // namespace

TEST_P(TermAccumulation, LikeRealCoefficientsKeepTheUnitRemainderInEveryOrder) {
    auto site = mdl::spin();
    std::array<double, 3> values{-A, 1., A};
    do {
        SCOPED_TRACE(::testing::PrintToString(values));
        mdl::TermSum terms;
        for (auto value : values) terms = terms + value * op(site, "Sz");
        expect_action(lower(GetParam(), terms, site, 1, {"Sz"}), site,
                      {.5, 0., 0., -.5});
    } while (std::next_permutation(values.begin(), values.end()));
}

TEST_P(TermAccumulation, LikeComplexCoefficientsKeepBothComponentsInEveryOrder) {
    auto site = mdl::spin();
    std::array<C, 3> const values{C(A, -A), C(.75, .25), C(-A, A)};
    std::array<int, 3> order{0, 1, 2};
    do {
        SCOPED_TRACE(::testing::PrintToString(order));
        mdl::TermSum terms;
        for (auto index : order) terms = terms + values[static_cast<std::size_t>(index)] * op(site, "Sz");
        // Generic operator compilation must not impose Hermiticity.
        expect_action(lower(GetParam(), terms, site, 1, {"Sz"}), site,
                      {C(.375, .125), 0., 0., C(-.375, -.125)});
    } while (std::next_permutation(order.begin(), order.end()));
}

TEST_P(TermAccumulation, ExactCancellationDoesNotLeaveAnOrderDependentOperator) {
    auto site = mdl::spin();
    std::array<double, 4> values{-A, -1., 1., A};
    do {
        SCOPED_TRACE(::testing::PrintToString(values));
        mdl::TermSum terms;
        for (auto value : values) terms = terms + value * op(site, "Sz");
        expect_action(lower(GetParam(), terms, site, 1, {"Sz"}), site, std::vector<C>(4));
    } while (std::next_permutation(values.begin(), values.end()));
}

TEST_P(TermAccumulation, FiniteScalarSurvivesOverflowingIntermediateCoefficientSums) {
    auto site = mdl::spin();
    for (auto const& values : std::vector<std::vector<double>>{
             {B, B, -B, -B, 1.}, {B, 1., B, -B, -B}, {-B, -B, 1., B, B}}) {
        SCOPED_TRACE(::testing::PrintToString(values));
        mdl::TermSum terms;
        for (auto value : values) terms = terms + value * op(site, "I");
        expect_action(lower(GetParam(), terms, site), site, {1., 0., 0., 1.});
    }
}

TEST_P(TermAccumulation, UnrepresentableMergedCoefficientDoesNotRejectFiniteSpinMatrix) {
    auto site = mdl::spin();
    auto terms = B * op(site, "Sz") + B * op(site, "Sz");
    // 2B cannot be a double coefficient, but (2B) Sz has entries +/-B.
    // Keeping finite contributions or moving scale into factors are both valid.
    expect_action(lower(GetParam(), terms, site, 1, {"Sz"}), site, {B, 0., 0., -B});
}

TEST_P(TermAccumulation, OrderedOddFactorsRetainGradedSignsAndComplexRemainders) {
    auto site = mdl::spinless_fermion();
    auto forward = op(site, "c_dag", 2) * op(site, "c", 0);
    auto reversed = op(site, "c", 0) * op(site, "c_dag", 2);
    C const remainder(.75, .25);
    for (auto convention : {strict, mdl::FermionConvention::lattice}) {
        // The QMC descriptor has its own JW wire convention even when
        // algebraic swap signs are disabled. Its lattice interpretation is
        // outside this arithmetic correction; lock strict QMC and both
        // conventions of the general Hamiltonian compiler here.
        if (GetParam() == Route::local && convention == mdl::FermionConvention::lattice) continue;
        SCOPED_TRACE(convention == strict ? "strict" : "lattice");
        // Reversing distinct odd factors gives -T under CAR, +T under
        // the lattice tensor convention. The complete operator is rT+h.c.
        std::array<mdl::TermSum, 3> const pieces{
            A * forward, remainder * forward, (convention == strict ? A : -A) * reversed};
        std::array<int, 3> order{0, 1, 2};
        do {
            SCOPED_TRACE(::testing::PrintToString(order));
            mdl::TermSum terms;
            for (auto index : order) terms = terms + pieces[static_cast<std::size_t>(index)];
            terms = terms.add_hc();
            EXPECT_NO_THROW(mdl::validate_hermiticity(terms, {site, site, site}, convention, 0.));
            std::vector<C> expected(64);
            expected[1 * 8 + 4] = remainder;
            expected[4 * 8 + 1] = std::conj(remainder);
            C const crossed = convention == strict ? -remainder : remainder;
            expected[3 * 8 + 6] = crossed;
            expected[6 * 8 + 3] = std::conj(crossed);
            expect_action(lower(GetParam(), terms, site, 3, {"N"}, convention), site, expected);
        } while (std::next_permutation(order.begin(), order.end()));
    }
}

TEST_P(TermAccumulation, SameSiteProductsRetainElementaryTransitionAndOrder) {
    for (C remainder : {C(1.), C(.75, -.25)}) {
        SCOPED_TRACE(remainder);
        auto site = product_site(remainder);
        std::vector<C> expected(9);
        expected[1] = remainder;
        expect_action(lower(GetParam(), op(site, "X") * op(site, "Y"), site), site, expected);
        // YX=0, independently of the cancellation in XY. Product order matters.
        expect_action(lower(GetParam(), op(site, "Y") * op(site, "X"), site), site, std::vector<C>(9));
    }
}

TEST_P(TermAccumulation, SameSiteProductIsValidatedBeforeFalseChargeProjection) {
    auto site = product_site();
    auto terms = op(site, "X") * op(site, "Y");
    auto before = terms;
    EXPECT_THROW(lower(GetParam(), terms, site, 1, {"Q"}), std::invalid_argument);
    EXPECT_THROW(mdl::compile_local_terms(terms, {site}, {"Q"}), std::invalid_argument);
    expect_unchanged(terms, before);
}

TEST_P(TermAccumulation, HermitianCombinationOfSameSiteProductsRetainsBothTransitions) {
    C const remainder(.75, -.25);
    auto site = product_site(remainder);
    auto terms = (op(site, "X") * op(site, "Y")).add_hc();
    EXPECT_NO_THROW(mdl::validate_hermiticity(terms, {site}, strict, 0.));
    std::vector<C> expected(9);
    expected[1] = remainder;
    expected[3] = std::conj(remainder);
    expect_action(lower(GetParam(), terms, site), site, expected);
    EXPECT_THROW(lower(GetParam(), terms, site, 1, {"Q"}), std::invalid_argument);
}

TEST_P(TermAccumulation, ThreeFactorProductKeepsIntermediateCancellationUntilFinalProjector) {
    auto site = product_site();
    auto terms = op(site, "X") * op(site, "Y") * op(site, "Z");
    // XY=E01 and Z=E10, hence XYZ=E00: Hermitian and Q-conserving.
    EXPECT_NO_THROW(mdl::validate_hermiticity(terms, {site}, strict, 0.));
    std::vector<C> expected(9);
    expected[0] = 1.;
    expect_action(lower(GetParam(), terms, site, 1, {"Q"}), site, expected);
}

TEST_P(TermAccumulation, ThreeFactorProductCanReturnFromIntermediateOverflowToFiniteMatrix) {
    mdl::SiteBasis description;
    description.dim = 2;
    description.quantum_numbers["Q"] = {0., 1.};
    description.operators["X"] = {"X", 2, {B, 0., 0., 0.}};
    description.operators["Y"] = {"Y", 2, {2., 0., 0., 1.}};
    description.operators["Z"] = {"Z", 2, {.5, 0., 0., 1.}};
    auto site = mdl::finalize_site(std::move(description));
    auto terms = op(site, "X") * op(site, "Y") * op(site, "Z");
    // XY contains 2B, but XYZ is exactly B E00. Rounding/materializing
    // the intermediate matrix is not a valid reason to reject this input.
    EXPECT_NO_THROW(mdl::validate_hermiticity(terms, {site}, strict, 0.));
    expect_action(lower(GetParam(), terms, site, 1, {"Q"}), site, {B, 0., 0., 0.});
}

TEST_P(TermAccumulation, ExternalCoefficientCompensatesOverflowAndUnderflowOfLocalProduct) {
    struct Case { double factor, coefficient, expected; };
    for (auto const& sample : {Case{0x1p700, 0x1p-400, 0x1p1000},
                               Case{0x1p-700, 0x1p400, 0x1p-1000}}) {
        SCOPED_TRACE(sample.factor);
        mdl::SiteBasis description;
        description.dim = 2;
        description.quantum_numbers["Q"] = {0., 1.};
        description.operators["P"] = {"P", 2, {sample.factor, 0., 0., 0.}};
        auto site = mdl::finalize_site(std::move(description));
        auto terms = sample.coefficient * (op(site, "P") * op(site, "P"));
        // P^2 alone is outside double range (or below its smallest subnormal),
        // but the supplied complete term is an exactly representable projector.
        EXPECT_NO_THROW(mdl::validate_hermiticity(terms, {site}, strict, 0.));
        expect_action(lower(GetParam(), terms, site, 1, {"Q"}), site,
                      {sample.expected, 0., 0., 0.});
    }
}

TEST_P(TermAccumulation, PartiallyRoundedSubnormalProductCannotBeAmplifiedAsExactInput) {
    mdl::SiteBasis description;
    description.dim = 2;
    description.quantum_numbers["Q"] = {0., 1.};
    description.operators["X"] = {"X", 2, {0x1.8p-537, 0., 0., 0.}};
    description.operators["Y"] = {"Y", 2, {0x1p-537, 0., 0., 0.}};
    auto site = mdl::finalize_site(std::move(description));
    auto terms = 0x1p10 * (op(site, "X") * op(site, "Y"));
    // (3*2^-538)*(2^-537)*(2^10) = 3*2^-1065 exactly. The unscaled
    // product rounds 1.5 subnormal units to 2; it is nonzero but inaccurate.
    expect_action(lower(GetParam(), terms, site, 1, {"Q"}), site,
                  {0x1.8p-1064, 0., 0., 0.});
}

TEST_P(TermAccumulation, ZeroEliminationDoesNotHideMalformedOrNonfiniteFactors) {
    auto site = mdl::spin();
    mdl::Term bad_shape;
    bad_shape.coeff = 0.;
    bad_shape.factors.push_back({0, {"bad_shape", 2, {1., 0., 0.}}, site});
    mdl::TermSum const shape_terms({bad_shape});
    EXPECT_THROW(lower(GetParam(), shape_terms, site), std::invalid_argument);

    mdl::Term bad_value;
    bad_value.factors.push_back({0, {"zero", 2, {0., 0., 0., 0.}}, site});
    bad_value.factors.push_back({0, {"nonfinite", 2,
        {std::numeric_limits<double>::infinity(), 0., 0., 0.}}, site});
    mdl::TermSum const value_terms({bad_value});
    EXPECT_THROW(lower(GetParam(), value_terms, site), std::invalid_argument);
    expect_unchanged(shape_terms, mdl::TermSum({bad_shape}));
    expect_unchanged(value_terms, mdl::TermSum({bad_value}));
}

TEST_P(TermAccumulation, CompilationPreservesSourcesAndOwnsTheResult) {
    auto writable = std::make_shared<mdl::SiteBasis>(*product_site());
    mdl::SiteBasisPtr site = writable;
    auto terms = op(site, "X") * op(site, "Y") * op(site, "Z");
    auto before = terms;
    auto const original_x = site->op("X").matrix;
    auto const input = lower(GetParam(), terms, site, 1, {"Q"});
    expect_unchanged(terms, before);
    EXPECT_EQ(site->op("X").matrix, original_x);
    EXPECT_EQ(input.diag_observables.at(0).at("Q"), (std::vector<double>{0., 1., 2.}));
    terms = mdl::TermSum::zero();
    before = mdl::TermSum::zero();
    writable->operators["X"].matrix.assign(9, 17.);
    writable->quantum_numbers["Q"][0] = 19.;
    std::vector<C> expected(9);
    expected[0] = 1.;
    expect_action(input, site, expected);
    EXPECT_EQ(input.diag_observables.at(0).at("Q"), (std::vector<double>{0., 1., 2.}));
}

TEST(TermAccumulationValidation, NonHermitianSameSiteProductCannotDisappearBeforeValidation) {
    for (C remainder : {C(1.), C(.75, -.25)}) {
        auto site = product_site(remainder);
        auto terms = op(site, "X") * op(site, "Y");
        auto before = terms;
        EXPECT_THROW(mdl::validate_hermiticity(terms, {site}), std::invalid_argument);
        expect_unchanged(terms, before);
    }
}

INSTANTIATE_TEST_SUITE_P(PublicCompilers, TermAccumulation,
    ::testing::Values(Route::hamiltonian, Route::local),
    [](auto const& info) { return info.param == Route::hamiltonian ? "Hamiltonian" : "Local"; });
