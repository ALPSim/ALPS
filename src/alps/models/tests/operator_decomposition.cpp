// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/operator_decomposition.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// dense comparisons use the LAPACK-based test helper instead of Eigen.
#include <alps/models/hamiltonian.hpp>
#include <gtest/gtest.h>
#include "dense_e0.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace mdl = alps::models;

TEST(OperatorDecomposition, KeepsChargeChangingHamiltonianAndObservableTerms) {
    auto site = mdl::spin();
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    auto flip = mdl::TermSum::located(site, "Sx", 0);
    h.add(flip);
    // Deliberate API correction: MC lowering rejects a false declaration;
    // the unrestricted operator/observable route still retains every term.
    EXPECT_THROW(h.to_local_decomposition(), std::invalid_argument);
    auto input = h.operator_decomposition();
    ASSERT_EQ(input.op.terms.size(), 1);
    EXPECT_EQ(input.op.terms[0].factors[0].matrix, site->op("Sx").matrix);
    auto observable = h.operator_decomposition(mdl::TermSum::located(site, "Sp", 1));
    ASSERT_EQ(observable.op.terms.size(), 1);
    ASSERT_EQ(observable.op.terms[0].factors.size(), 1);
    EXPECT_EQ(observable.op.terms[0].factors[0].site, 1);
    EXPECT_EQ(observable.op.terms[0].factors[0].matrix, site->op("Sp").matrix);
    EXPECT_EQ(observable.op.site_dims, (std::vector<int>{2, 2}));
    EXPECT_EQ(observable.diag_observables.at(0).at("Sz"), site->quantum_numbers.at("Sz"));
}

TEST(OperatorDecomposition, StrictAndLatticeFermionConventions) {
    auto site = mdl::spinless_fermion();
    auto op = mdl::TermSum::located(site, "c", 1) * mdl::TermSum::located(site, "c", 0);
    mdl::Hamiltonian strict(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"N"});
    mdl::Hamiltonian local(alps::lattice::chain(2), site, mdl::FermionConvention::lattice, {"N"});
    auto a = strict.operator_decomposition(op);
    auto b = local.operator_decomposition(op);
    ASSERT_EQ(a.op.terms.size(), 1);
    ASSERT_EQ(b.op.terms.size(), 1);
    EXPECT_EQ(a.op.terms[0].coefficient, mdl::Complex(-1.));
    EXPECT_EQ(b.op.terms[0].coefficient, mdl::Complex(1.));
    EXPECT_EQ(a.op.terms[0].jordan_wigner, (std::vector<int>{0}));
    EXPECT_TRUE(b.op.terms[0].jordan_wigner.empty());
}

TEST(OperatorDecompositionConservation, RejectsDeclaredChargeBrokenByAnOnsiteField) {
    auto site = mdl::spin();
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    h.add(mdl::TermSum::located(site, "Sx", 0));
    // A declared symmetry is a claim to validate, not an implicit projection.
    EXPECT_THROW(h.to_local_decomposition(), std::invalid_argument);
    auto complete = alps_test::dense_hamiltonian(h, h.operator_decomposition());
    EXPECT_EQ(complete(0, 2), mdl::Complex(.5));
    EXPECT_EQ(complete(2, 0), mdl::Complex(.5));
}

TEST(OperatorDecompositionConservation, RejectsUncancelledDoubleSpinFlips) {
    auto site = mdl::spin();
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    h.add(op("Sx", 0) * op("Sx", 1) + .75 * op("Sy", 0) * op("Sy", 1));
    // The all-up/all-down coupling is (1 - .75)/4, not a roundoff residual.
    auto complete = alps_test::dense_hamiltonian(h, h.operator_decomposition());
    EXPECT_EQ(complete(0, 3), mdl::Complex(1. / 16.));
    EXPECT_THROW(h.to_local_decomposition(), std::invalid_argument);
}

TEST(OperatorDecompositionConservation, CancellingCartesianSpinTermsRemainValid) {
    auto site = mdl::spin();
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    h.add(op("Sx", 0) * op("Sx", 1) + op("Sy", 0) * op("Sy", 1));
    auto expected = alps_test::DenseMatrix(4);
    expected(1, 2) = expected(2, 1) = .5;
    auto actual = alps_test::dense_hamiltonian(h, h.to_local_decomposition());
    EXPECT_EQ(alps_test::frobenius_distance(actual, expected), 0.);
}

TEST(OperatorDecompositionConservation, ComplexCoefficientsRetainConservationAndHermitianAction) {
    auto site = mdl::spin();
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    mdl::Complex const coefficient(.75, -.25);
    h.add(coefficient * op("Sp", 0) * op("Sm", 1), true);
    auto expected = alps_test::DenseMatrix(4);
    expected(1, 2) = coefficient;
    expected(2, 1) = std::conj(coefficient);
    auto actual = alps_test::dense_hamiltonian(h, h.to_local_decomposition());
    EXPECT_EQ(alps_test::frobenius_distance(actual, expected), 0.);

    mdl::Hamiltonian broken(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    broken.add(coefficient * op("Sp", 0), true);
    EXPECT_THROW(broken.to_local_decomposition(), std::invalid_argument);
}

TEST(OperatorDecompositionConservation, ChargeChangingTermsCancelAcrossDifferentSupports) {
    mdl::SiteBasis basis = *mdl::spin();
    basis.operators["Pup"] = {"Pup", 2, {1., 0., 0., 0.}};
    basis.operators["Pdown"] = {"Pdown", 2, {0., 0., 0., 1.}};
    auto site = mdl::finalize_site(std::move(basis));
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    h.add(op("Sx", 0) * op("Pup", 1) + op("Sx", 0) * op("Pdown", 1) - op("Sx", 0));
    h.add(op("Sz", 0) * op("Sz", 1));
    auto expected = alps_test::DenseMatrix(4);
    expected(0, 0) = .25;
    expected(1, 1) = -.25;
    expected(2, 2) = -.25;
    expected(3, 3) = .25;
    auto actual = alps_test::dense_hamiltonian(h, h.to_local_decomposition());
    EXPECT_EQ(alps_test::frobenius_distance(actual, expected), 0.);
}

TEST(OperatorDecompositionConservation, ValidationDoesNotConstructTheGlobalHilbertSpace) {
    auto site = mdl::spin();
    mdl::Hamiltonian h(alps::lattice::chain(64), site, mdl::FermionConvention::strict, {"Sz"});
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    h.add_bond(1., [&](int i, int j) {
        return op("Sx", i) * op("Sx", j) + op("Sy", i) * op("Sy", j);
    });
    // A full 2^64-dimensional operator cannot be built for validation.
    auto local = h.to_local_decomposition();
    EXPECT_EQ(local.op.site_dims.size(), 64);
    EXPECT_FALSE(local.op.terms.empty());
    for (auto const& term : local.op.terms) EXPECT_LE(term.factors.size(), 2);
}

TEST(OperatorDecompositionAmplitude, QuantumNumberSplittingPreservesTinyComplexEntriesExactly) {
    mdl::LocalOp op{"mixed", 3,
        {2., mdl::Complex(1e-12, -2e-12), 0.,
         mdl::Complex(-3e-13, 4e-13), 0., 1.,
         0., 1., -1.}};
    std::map<std::string, std::vector<double>> qns{{"N", {0., 1., 2.}}};
    auto parts = op.split_by_qn_delta(qns);
    mdl::Matrix reconstructed(op.matrix.size());
    for (auto const& part : parts)
        for (std::size_t i = 0; i < reconstructed.size(); ++i)
            reconstructed[i] += part.matrix[i];
    // Splitting merely partitions entries. No numerical approximation is
    // needed, including for entries below the charge-label tolerance.
    EXPECT_EQ(reconstructed, op.matrix);
}

TEST(OperatorDecompositionAmplitude, ChargeInferenceIncludesAmplifiableSmallEntries) {
    mdl::Matrix mixed{0., 1., mdl::Complex(1e-12, -2e-12), 0.};
    auto deltas = mdl::compute_qn_deltas(mixed, 2, {{"Sz", {.5, -.5}}});
    EXPECT_FALSE(deltas.count("Sz") != 0);
    mdl::Matrix raising{0., mdl::Complex(1e-12, -2e-12), 0., 0.};
    auto raising_deltas = mdl::compute_qn_deltas(raising, 2, {{"Sz", {.5, -.5}}});
    ASSERT_TRUE(raising_deltas.count("Sz") != 0);
    EXPECT_EQ(raising_deltas.at("Sz"), 1.);
}

TEST(OperatorDecompositionAmplitude, ChargeInferenceOmitsNonComparableDeltas) {
    mdl::Matrix raising{0., 1., 0., 0.};
    double const nan = std::numeric_limits<double>::quiet_NaN();
    double const inf = std::numeric_limits<double>::infinity();
    for (double charge : {nan, inf})
        EXPECT_TRUE(mdl::compute_qn_deltas(raising, 2, {{"Q", {charge, 0.}}}).empty());
    for (double tolerance : {-1., nan})
        EXPECT_TRUE(mdl::compute_qn_deltas(raising, 2, {{"Q", {0., 1.}}}, tolerance).empty());
    // No transitions still preserve every declared charge, independently of tolerance.
    EXPECT_EQ(mdl::compute_qn_deltas(mdl::Matrix(4), 2, {{"Q", {nan, inf}}}, nan),
              (std::map<std::string, double>{{"Q", 0.}}));
}

TEST(OperatorDecompositionAmplitude, MovingScaleBetweenCoefficientAndLocalOperatorPreservesAction) {
    // One mixed-charge operator has a large off-diagonal component and a
    // tiny diagonal component. The off-diagonal parts cancel in X + X^dag,
    // leaving an order-one conserving diagonal after coefficient scaling.
    mdl::SiteBasis basis;
    basis.dim = 2;
    basis.quantum_numbers["N"] = {0., 1.};
    basis.operators["X"] = {"X", 2,
        {mdl::Complex(1e-12, -2e-12), 1., -1., mdl::Complex(-2e-12, 3e-12)}};
    auto small_site = mdl::finalize_site(basis);
    for (auto& entry : basis.operators["X"].matrix) entry *= 1e12;
    auto scaled_site = mdl::finalize_site(std::move(basis));
    auto expected = alps_test::DenseMatrix(4);
    expected(0, 0) = 2.;
    expected(1, 1) = 2.;
    expected(2, 2) = -4.;
    expected(3, 3) = -4.;
    for (auto const& [site, coefficient] :
         std::vector<std::pair<mdl::SiteBasisPtr, double>>{{small_site, 1e12}, {scaled_site, 1.}}) {
        mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"N"});
        h.add(coefficient * mdl::TermSum::located(site, "X", 0), true);
        auto actual = alps_test::dense_hamiltonian(h, h.to_local_decomposition());
        EXPECT_NEAR(alps_test::frobenius_distance(actual, expected), 0., 1e-14);
    }
}

TEST(OperatorDecompositionAmplitude, NearIdentityFactorsKeepAmplifiedOffDiagonalEntries) {
    mdl::SiteBasis basis = *mdl::spin();
    basis.operators["NearI"] = {"NearI", 2, {1., 1e-14, 1e-14, 1.}};
    auto site = mdl::finalize_site(std::move(basis));
    auto term = 1e14 * mdl::TermSum::located(site, "NearI", 0);
    mdl::Hamiltonian h(alps::lattice::chain(2), site);
    h.add(term);
    auto actual = alps_test::dense_hamiltonian(h, h.to_local_decomposition());
    EXPECT_NEAR(std::abs(actual(0, 2) - mdl::Complex(1.)), 0., 1e-14);
    EXPECT_NEAR(std::abs(actual(2, 0) - mdl::Complex(1.)), 0., 1e-14);
    EXPECT_EQ(actual(0, 0), mdl::Complex(1e14));

    mdl::Hamiltonian broken(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    broken.add(term);
    EXPECT_THROW(broken.to_local_decomposition(), std::invalid_argument);
}

TEST(OperatorDecompositionConservation, LongRangeFermionChargeValidationKeepsLocalCost) {
    auto site = mdl::spinless_fermion();
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    for (auto convention : {mdl::FermionConvention::strict, mdl::FermionConvention::lattice}) {
        mdl::Hamiltonian hopping(alps::lattice::chain(64), site, convention, {"N"});
        hopping.add(mdl::Complex(.75, .25) * op("c_dag", 63) * op("c", 0), true);
        EXPECT_FALSE(hopping.to_local_decomposition().op.terms.empty());

        mdl::Hamiltonian pairing(alps::lattice::chain(64), site, convention, {"N"});
        pairing.add(op("c", 63) * op("c", 0), true);
        EXPECT_THROW(pairing.to_local_decomposition(), std::invalid_argument);
    }
}

TEST(OperatorDecompositionConservation, LargeCancellingTermsCannotHideARepresentableRemainder) {
    mdl::SiteBasis basis = *mdl::spin();
    basis.operators["Pup"] = {"Pup", 2, {1., 0., 0., 0.}};
    basis.operators["Pdown"] = {"Pdown", 2, {0., 0., 0., 1.}};
    auto site = mdl::finalize_site(std::move(basis));
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    double const amplitude = std::ldexp(1., 50);
    for (double remainder : {-1., 0., 1.}) for (bool reverse : {false, true}) {
        SCOPED_TRACE(remainder);
        SCOPED_TRACE(reverse);
        std::vector<mdl::TermSum> pieces{
            amplitude * op("Sx", 0) * op("Pup", 1),
            amplitude * op("Sx", 0) * op("Pdown", 1),
            -(amplitude - remainder) * op("Sx", 0)};
        if (reverse) std::reverse(pieces.begin(), pieces.end());
        mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
        for (auto const& piece : pieces) h.add(piece);
        // Pup + Pdown = I, so the exact operator is remainder*Sx(0).
        // Its off-diagonal entries are exactly +/- .5, regardless of the
        // much larger coefficients that cancelled to produce it.
        if (remainder != 0.) EXPECT_THROW(h.to_local_decomposition(), std::invalid_argument);
        else EXPECT_NO_THROW(h.to_local_decomposition());
    }
}

TEST(OperatorDecompositionConservation, ValidateBeforeDoubleLikeTermAccumulation) {
    auto site = mdl::spin();
    double const amplitude = std::ldexp(1., 54);
    for (auto const& coefficients : std::vector<std::vector<double>>{
             {amplitude, 1., -amplitude}, {amplitude, -amplitude, 1.}}) {
        mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
        for (double coefficient : coefficients)
            h.add(coefficient * mdl::TermSum::located(site, "Sx", 0));
        // The supplied symbolic sum is Sx(0); ordinary double left-to-right
        // coefficient accumulation loses its unit coefficient in one order.
        EXPECT_THROW(h.to_local_decomposition(), std::invalid_argument);
    }
}

TEST(OperatorDecompositionConservation, PerFactorChargeToleranceDoesNotAccumulateSilently) {
    mdl::SiteBasis a;
    a.dim = 3;
    a.quantum_numbers["Q"] = {0., 1., 2. + .75e-9};
    a.operators["R"] = {"R", 3, {0., 0., 0., 1., 0., 0., 0., 1., 0.}};
    mdl::SiteBasis b;
    b.dim = 3;
    b.quantum_numbers["Q"] = {0., 1., 2. - .75e-9};
    b.operators["L"] = {"L", 3, {0., 1., 0., 0., 0., 1., 0., 0., 0.}};
    auto first = mdl::finalize_site(std::move(a));
    auto second = mdl::finalize_site(std::move(b));
    auto term = mdl::TermSum::located(first, "R", 0) * mdl::TermSum::located(second, "L", 1);
    // The E21 tensor E12 entry changes total Q by 1.5e-9, although each
    // factor's legacy inferred label lies within .75e-9 of +/-1.
    EXPECT_THROW(mdl::compile_local_terms(term, {first, second}, {"Q"}), std::invalid_argument);
}

TEST(OperatorDecompositionConservation, SplittingCannotReintroduceCancelledTransitionsNearTolerance) {
    mdl::SiteBasis basis;
    basis.dim = 3;
    basis.quantum_numbers["Q"] = {0., .75e-9, 2.25e-9};
    basis.operators["R"] = {"R", 3, {0., 0., 0., 1., 0., 0., 0., 1., 0.}};
    basis.operators["B"] = {"B", 3, {0., 0., 0., 0., 0., 0., 0., 1., 0.}};
    auto site = mdl::finalize_site(std::move(basis));
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Q"});
    h.add(mdl::TermSum::located(site, "R", 0) - mdl::TermSum::located(site, "B", 0), true);
    auto expected = alps_test::DenseMatrix(9);
    for (int spectator = 0; spectator < 3; ++spectator) {
        expected(3 + spectator, spectator) = 1.;
        expected(spectator, 3 + spectator) = 1.;
    }
    auto actual = alps_test::dense_hamiltonian(h, h.to_local_decomposition());
    EXPECT_EQ(alps_test::frobenius_distance(actual, expected), 0.);
}

TEST(OperatorDecompositionConservation, RoundedComponentLabelsCannotChangeTheKeptOperator) {
    mdl::SiteBasis basis;
    basis.dim = 2;
    basis.quantum_numbers["Q"] = {0., 1.0004e-9};
    basis.operators["R"] = {"R", 2, {1., 0., 1., 0.}};
    basis.operators["B"] = {"B", 2, {0., 0., 1., 0.}};
    auto site = mdl::finalize_site(std::move(basis));
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Q"});
    h.add(mdl::TermSum::located(site, "R", 0) - mdl::TermSum::located(site, "B", 0));
    auto expected = alps_test::DenseMatrix(4);
    expected(0, 0) = expected(1, 1) = 1.;
    auto actual = alps_test::dense_hamiltonian(h, h.to_local_decomposition());
    EXPECT_EQ(alps_test::frobenius_distance(actual, expected), 0.);
}

TEST(OperatorDecompositionConservation, MissingLocalChargesMeanNeutralNotMissingGlobalSymmetry) {
    auto fermion = mdl::spinless_fermion();
    auto spin = mdl::spin();
    auto neutral = mdl::TermSum::located(fermion, "n", 0) * mdl::TermSum::located(spin, "Sx", 1);
    EXPECT_NO_THROW(mdl::compile_local_terms(neutral, {fermion, spin}, {"N"}));
    auto charged = mdl::TermSum::located(fermion, "c", 0) * mdl::TermSum::located(spin, "Sx", 1);
    EXPECT_THROW(mdl::compile_local_terms(charged, {fermion, spin}, {"N"}), std::invalid_argument);
    EXPECT_THROW(mdl::compile_local_terms(neutral, {fermion, spin}, {"Unknown"}), std::invalid_argument);
}
