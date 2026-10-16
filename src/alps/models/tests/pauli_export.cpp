// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Sparse Pauli export of spin and boson models against angular-momentum
// matrices written independently of the kit, the trace-completed local
// identity on unused code words, projectors and penalties, pruning bounds,
// resource limits and stale charge declarations. Cases follow
// skilledwolf/alpsx (branch agent/v3-execution-and-api-audit at f8ee6e98b,
// src/model/tests/models/pauli_export.cpp and tests/alps/test_spin_pauli.py),
// adapted to GoogleTest, the CompiledOperator/PhysicalBasis entry and LAPACK
// eigenvalues (alpsx used Eigen). The dense Kronecker assembly of
// dense_e0.hpp supplies an independent spectrum.
#include "pauli_support.hpp"
#include <alps/models/builtin.hpp>
#include <alps/models/hamiltonian.hpp>
#include <alps/models/pauli.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;
using namespace alps_test;
using mdl::QubitEncoding;
using mdl::QubitEncodingScheme;

namespace {
constexpr QubitEncodingScheme schemes[] = {QubitEncodingScheme::binary, QubitEncodingScheme::gray};

mdl::TermSum L(mdl::SiteBasisPtr const& site, char const* name, int i) { return mdl::TermSum::located(site, name, i); }

// Arbitrary observables: no Hermiticity, parity or charge requirement, no pruning.
mdl::PauliExportOptions observe() {
    mdl::PauliExportOptions options;
    options.require_hermitian = options.require_even = options.require_conservation = false;
    options.tolerance = 0.;
    return options;
}
mdl::PauliOperator observable(mdl::Hamiltonian const& model, mdl::TermSum const& op, QubitEncoding const& encoding) {
    return mdl::export_pauli(model.operator_decomposition(op).op, mdl::physical_basis(model), encoding, observe());
}
mdl::PauliExportOptions exact_options() {
    mdl::PauliExportOptions options;
    options.tolerance = 0.;
    options.validation_tolerance = 1e-12;
    return options;
}
template <class F>
void expect_invalid(F&& f, std::string const& fragment) {
    try {
        f();
        FAIL() << "expected invalid_argument containing: " << fragment;
    } catch (std::invalid_argument const& e) {
        EXPECT_NE(std::string(e.what()).find(fragment), std::string::npos) << e.what();
    }
}
DenseMatrix kron3(DenseMatrix const& a, DenseMatrix const& b, DenseMatrix const& c) { return kron(kron(a, b), c); }
} // namespace

TEST(PauliExport, SpinHalfDimerWithComplexOperatorAndCasimir) {
    auto const spin = mdl::spin();
    mdl::Hamiltonian model(lat::chain(2, 1.0, lat::Boundary::open), spin);
    model.add(.7 * L(spin, "Sy", 0) * L(spin, "Sx", 1) + .13 * L(spin, "Sz", 1) + .19 * L(spin, "I", 0));
    auto const [x, y, z] = spin_matrices(.5);
    auto const i2 = DenseMatrix::identity(2);
    auto expected = add(add(scale(kron(y, x), .7), kron(i2, z), .13), DenseMatrix::identity(4), .19);
    auto const exported = mdl::export_pauli(model);
    EXPECT_EQ(exported.num_qubits, 2);
    auto const encoding = mdl::make_qubit_encoding(model);
    // Two-state sites fill their qubit, but site 0 is the least significant
    // one: the physical index map, not the raw register order, is the contract.
    auto const register_matrix = pauli_matrix(exported);
    EXPECT_LT(max_abs_difference(restrict_to(register_matrix, physical_indices(encoding)), expected), 1e-14);
    EXPECT_GT(max_abs_difference(register_matrix, expected), .1);
    auto const casimir = L(spin, "Sx", 0) * L(spin, "Sx", 0) + L(spin, "Sy", 0) * L(spin, "Sy", 0)
                       + L(spin, "Sz", 0) * L(spin, "Sz", 0);
    EXPECT_LT(max_abs_difference(pauli_matrix(observable(model, casimir, encoding)),
                                 scale(DenseMatrix::identity(4), .75)), 1e-14);
}

TEST(PauliExport, EncodedDimerActionObservablesAndStates) {
    for (double S : {.5, 1., 1.5, 2.})
        for (auto scheme : schemes) {
            SCOPED_TRACE("S=" + std::to_string(S) + (scheme == QubitEncodingScheme::gray ? " gray" : " binary"));
            auto const model = mdl::heisenberg(lat::chain(2, 1.0, lat::Boundary::open), 1., 0., S);
            auto const site = model.site_for(0);
            auto const encoding = mdl::make_qubit_encoding(model, scheme);
            auto const exported = mdl::export_pauli(model, encoding, exact_options());
            auto const [x, y, z] = spin_matrices(S);
            auto const d = x.dim;
            auto const identity = DenseMatrix::identity(d);
            auto const reference = add(add(kron(x, x), kron(y, y)), kron(z, z));
            auto const indices = physical_indices(encoding);
            auto const actual = pauli_matrix(exported);
            EXPECT_LT(action_distance(actual, indices, reference), 1e-13);
            EXPECT_LT(max_abs_difference(actual, adjoint(actual)), 1e-14);
            std::mt19937 generator(823);
            std::normal_distribution<double> normal;
            std::vector<Complex> psi(d * d);
            double norm = 0.;
            for (auto& z_ : psi) { z_ = Complex(normal(generator), normal(generator)); norm += std::norm(z_); }
            for (auto& z_ : psi) z_ /= std::sqrt(norm);
            auto const mapped = encoding.encode_statevector(psi);
            for (std::size_t k = 0; k < psi.size(); ++k) EXPECT_EQ(mapped[indices[k]], psi[k]);
            auto const decoded = encoding.decode_statevector(mapped);
            EXPECT_EQ(decoded.amplitudes, psi);
            EXPECT_EQ(decoded.discarded_norm_squared, 0.);
            EXPECT_NEAR(encoding.physical_weight(mapped), 1., 1e-14);
            struct Case { mdl::TermSum op; DenseMatrix expected; };
            std::vector<Case> const cases{
                {L(site, "Sy", 0), kron(y, identity)},
                {L(site, "Sp", 1), kron(identity, add(x, y, Complex(0., 1.)))},
                {L(site, "Sz", 0) * L(site, "Sz", 0), kron(multiply(z, z), identity)}};
            for (auto const& [op, expected] : cases) {
                auto const matrix = pauli_matrix(observable(model, op, encoding));
                EXPECT_LT(action_distance(matrix, indices, expected), 1e-13);
                auto const lhs = vdot(mapped, matvec(matrix, mapped)), rhs = vdot(psi, matvec(expected, psi));
                EXPECT_LT(std::abs(lhs - rhs), 1e-13);
            }
        }
}

TEST(PauliExport, ArbitraryComplexMatrixAndIdentityCompletion) {
    for (int dimension : {3, 4, 5, 8, 9, 16})
        for (auto scheme : schemes) {
            SCOPED_TRACE(dimension);
            auto const site = mdl::spin((dimension - 1) / 2.);
            auto const basis = mdl::physical_basis({site});
            QubitEncoding const encoding(basis, scheme);
            std::mt19937 generator(static_cast<unsigned>(14 + dimension));
            std::normal_distribution<double> normal;
            auto const n = static_cast<std::size_t>(dimension);
            mdl::LocalFactor factor{0, std::vector<Complex>(n * n)};
            for (auto& z : factor.matrix) z = Complex(normal(generator), normal(generator));
            mdl::CompiledOperator const op{{dimension}, {{1., {factor}, {}}}};
            auto const actual = pauli_matrix(mdl::export_pauli(op, basis, encoding, observe()));
            auto const indices = physical_indices(encoding);
            Complex trace = 0.;
            for (std::size_t i = 0; i < n; ++i) trace += factor.matrix[i * n + i];
            // Unused code words carry the mean diagonal; the physical block is the matrix.
            DenseMatrix expected = scale(DenseMatrix::identity(actual.dim), trace / double(dimension));
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t c = 0; c < n; ++c)
                    expected(static_cast<std::size_t>(indices[r]), static_cast<std::size_t>(indices[c])) = factor.matrix[r * n + c];
            EXPECT_LT(max_abs_difference(actual, expected), 2e-14);
            mdl::LocalFactor identity{0, std::vector<Complex>(n * n)};
            for (std::size_t i = 0; i < n; ++i) identity.matrix[i * n + i] = 1.;
            auto const unit = mdl::export_pauli(mdl::CompiledOperator{{dimension}, {{1., {identity}, {}}}}, basis, encoding);
            EXPECT_EQ(unit.terms, (std::map<std::string, Complex>{{std::string(static_cast<std::size_t>(unit.num_qubits), 'I'), 1.}}));
            EXPECT_TRUE(mdl::export_pauli(mdl::CompiledOperator{{dimension}, {}}, basis, encoding).terms.empty());
        }
}

TEST(PauliExport, HeterogeneousBasisAndNoSpectatorProjectors) {
    std::map<int, mdl::SiteBasisPtr> const spaces{{0, mdl::spin(.5)}, {1, mdl::spin(1.)}, {2, mdl::spin(2.)}};
    lat::UnitCell cell(1);
    cell.add_vertex({0.}, 0);
    cell.add_vertex({.3}, 1);
    cell.add_vertex({.7}, 2);
    auto const lattice = cell.tile({1});
    ASSERT_EQ(lattice.num_sites(), 3);
    for (auto scheme : schemes) {
        mdl::Hamiltonian model(lattice, spaces);
        model.add(.31 * L(spaces.at(0), "Sy", 0) * L(spaces.at(2), "Sx", 2)
                  + .17 * L(spaces.at(1), "Sz", 1) * L(spaces.at(1), "Sz", 1));
        auto const encoding = mdl::make_qubit_encoding(model, scheme);
        EXPECT_EQ(encoding.num_qubits(), 6);
        auto const half = spin_matrices(.5), one = spin_matrices(1.), two = spin_matrices(2.);
        auto const expected = add(scale(kron3(half.y, DenseMatrix::identity(3), two.x), .31),
                                  kron3(DenseMatrix::identity(2), multiply(one.z, one.z), DenseMatrix::identity(5)), .17);
        EXPECT_LT(action_distance(pauli_matrix(mdl::export_pauli(model, encoding)), physical_indices(encoding), expected), 1e-13);
        mdl::Hamiltonian zero(lattice, spaces);
        zero.add(2.5 * L(spaces.at(0), "I", 0));
        EXPECT_EQ(mdl::export_pauli(zero, encoding).terms, (std::map<std::string, Complex>{{"IIIIII", 2.5}}));
    }
}

TEST(PauliExport, PhysicalBlockSpectrumMatchesTheDenseAssembly) {
    auto const model = mdl::heisenberg_d(lat::chain(3, 1.0, lat::Boundary::open), 1., .2, .1, 1.);
    auto const compiled = mdl::compile_hamiltonian(model);
    auto const reference = hermitian_eigenvalues(dense_hamiltonian(model, compiled));
    for (auto scheme : schemes) {
        auto const encoding = mdl::make_qubit_encoding(model, scheme);
        auto const exported = mdl::export_pauli(compiled.op, mdl::physical_basis(model), encoding, exact_options());
        auto const block = restrict_to(pauli_matrix(exported), physical_indices(encoding));
        auto const spectrum = hermitian_eigenvalues(block);
        ASSERT_EQ(spectrum.size(), reference.size());
        for (std::size_t k = 0; k < spectrum.size(); ++k) EXPECT_NEAR(spectrum[k], reference[k], 1e-12);
    }
}

TEST(PauliExport, ProjectorsPenaltyInvalidPatternsAndLostNorm) {
    for (auto scheme : schemes) {
        auto model = mdl::heisenberg(lat::chain(2, 1.0, lat::Boundary::open), 1., 0., 1.);
        auto const site = model.site_for(0);
        model.add(8. * L(site, "I", 0));
        auto const encoding = mdl::make_qubit_encoding(model, scheme);
        auto const encoded = mdl::export_pauli(model, encoding, exact_options());
        auto const indices = physical_indices(encoding);
        std::vector<std::uint64_t> invalid;
        for (std::uint64_t j = 0; j < 16; ++j)
            if (std::find(indices.begin(), indices.end(), j) == indices.end()) invalid.push_back(j);
        ASSERT_EQ(invalid.size(), 7u);
        for (auto j : invalid) {
            std::vector<int> bits;
            for (int b = 0; b < 4; ++b) bits.push_back(static_cast<int>((j >> b) & 1));
            EXPECT_FALSE(encoding.decode_basis_bits(bits).has_value());
        }
        auto const projector = multiply(pauli_matrix(mdl::local_valid_projector(encoding, 0)),
                                        pauli_matrix(mdl::local_valid_projector(encoding, 1)));
        DenseMatrix expected(16);
        for (auto j : indices) expected(static_cast<std::size_t>(j), static_cast<std::size_t>(j)) = 1.;
        EXPECT_EQ(max_abs_difference(projector, expected), 0.);
        auto const penalty = pauli_matrix(mdl::invalid_state_penalty(encoding, 2 * encoded.coefficient_l1(false) + 1));
        EXPECT_LT(max_column_norm(penalty, indices), 1e-12);
        auto const lifted = hermitian_eigenvalues(add(pauli_matrix(encoded), penalty));
        EXPECT_NEAR(lifted.front(), 6., 1e-12); // -2 + 8
        std::vector<Complex> state(16);
        state[indices[3]] = Complex(0., .6);
        state[invalid[0]] = .8;
        auto const decoded = encoding.decode_statevector(state);
        EXPECT_NEAR(decoded.discarded_norm_squared, .64, 1e-15);
        double kept = 0.;
        for (auto z : decoded.amplitudes) kept += std::norm(z);
        EXPECT_NEAR(kept, .36, 1e-15);
        EXPECT_NEAR(encoding.physical_weight(state), .36, 1e-15);
        EXPECT_THROW(mdl::invalid_state_penalty(encoding, -1.), std::invalid_argument);
        EXPECT_THROW(mdl::local_valid_projector(encoding, 2), std::invalid_argument);
    }
}

TEST(PauliExport, BasisSnapshotLimitsAndNoGlobalMatrix) {
    auto const model = mdl::heisenberg(lat::chain(2), 1., 0., 1.);
    auto const encoding = mdl::make_qubit_encoding(model);
    auto const compiled = mdl::compile_hamiltonian(model);
    auto reordered = mdl::physical_basis(model);
    std::reverse(reordered.sites[0].labels.begin(), reordered.sites[0].labels.end());
    expect_invalid([&] { (void)mdl::export_pauli(compiled.op, reordered, encoding); }, "basis mismatch");
    mdl::PauliExportOptions bounded;
    bounded.max_terms = 1;
    EXPECT_THROW(mdl::export_pauli(model, bounded), std::length_error);
    mdl::PauliExportOptions tiny;
    tiny.max_local_matrix_bytes = 63;
    EXPECT_THROW(mdl::export_pauli(model, tiny), std::length_error);
    EXPECT_THROW(mdl::export_pauli(model, encoding, tiny), std::length_error);
    // Just enough for a 4x4 complex workspace, too little for cached maps:
    // the uncached path, including charge validation, gives the same operator.
    mdl::PauliExportOptions uncached;
    uncached.max_local_matrix_bytes = 256;
    auto const reference = mdl::export_pauli(model);
    auto const sparse = mdl::export_pauli(model, uncached);
    EXPECT_EQ(sparse.terms, reference.terms);
    EXPECT_EQ(sparse.local_cache_hits, 0u);
    EXPECT_GT(reference.local_cache_hits, 0u);
    mdl::PauliExportOptions zero;
    zero.max_terms = 0;
    EXPECT_THROW(mdl::export_pauli(model, zero), std::invalid_argument);

    // A 128-site spin-2 chain decomposes three local matrices, never a global one.
    auto const two = mdl::spin(2.);
    mdl::Hamiltonian large(lat::chain(128, 1.0, lat::Boundary::open), two);
    large.add_bond(1., [two](int i, int j) { return L(two, "Sz", i) * L(two, "Sz", j); });
    auto const exported = mdl::export_pauli(large);
    EXPECT_EQ(exported.num_qubits, 384);
    EXPECT_GT(exported.local_cache_hits, 200u);
    EXPECT_LE(exported.local_decompositions, 3u);
    ASSERT_TRUE(exported.encoding.has_value());
    std::vector<int> const fours(128, 4);
    auto const bits = exported.encoding->encode_basis_bits(fours);
    ASSERT_EQ(bits.size(), 384u);
    auto const decoded = exported.encoding->decode_basis_bits(bits);
    ASSERT_TRUE(decoded.has_value());
    EXPECT_EQ(*decoded, fours);
    expect_invalid([&] { (void)exported.state_index(fours); }, "63 qubits");
}

TEST(PauliExport, PruningErrorBoundAndStaleCharge) {
    auto const site = mdl::spin(1.);
    mdl::Hamiltonian model(lat::chain(2), site);
    model.add(L(site, "Sx", 0) * L(site, "Sy", 1) + .003 * L(site, "Sz", 0));
    auto const exact = pauli_matrix(mdl::export_pauli(model, exact_options()));
    mdl::PauliExportOptions coarse;
    coarse.tolerance = .002;
    auto const approximate = mdl::export_pauli(model, coarse);
    EXPECT_GT(approximate.discarded_l1, 0.);
    auto const difference = hermitian_eigenvalues(add(pauli_matrix(approximate), exact, -1.));
    double const norm = std::max(std::abs(difference.front()), std::abs(difference.back()));
    EXPECT_GT(norm, 0.);
    EXPECT_LE(norm, approximate.discarded_l1 + 1e-14);
    // A stale conservation declaration is caught on the unpruned operator.
    auto const compiled = mdl::compile_hamiltonian(model);
    auto stale = mdl::physical_basis(model);
    stale.conserved_quantum_numbers = {"Sz"};
    expect_invalid([&] { (void)mdl::export_pauli(compiled.op, stale); }, "does not conserve declared quantum number Sz");
    mdl::PauliExportOptions pruned;
    pruned.tolerance = .9;
    pruned.validation_tolerance = 1e-12;
    expect_invalid([&] { (void)mdl::export_pauli(compiled.op, stale, pruned); }, "does not conserve");
    for (double bad : {-1., std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        mdl::PauliExportOptions options;
        options.validation_tolerance = bad;
        expect_invalid([&] { (void)mdl::export_pauli(model, options); }, "tolerance");
        options = {};
        options.tolerance = bad;
        expect_invalid([&] { (void)mdl::export_pauli(model, options); }, "tolerance");
    }
}

TEST(PauliExport, ExplicitBosonCutoffUsesTheSameFiniteSpaceContract) {
    auto const model = mdl::bose_hubbard(lat::chain(2, 1.0, lat::Boundary::open), .7, 1.3, 0., 0., 2);
    DenseMatrix a(3), n(3);
    a(0, 1) = 1.;
    a(1, 2) = std::sqrt(2.);
    n(1, 1) = 1.;
    n(2, 2) = 2.;
    auto const onsite = scale(multiply(n, add(n, DenseMatrix::identity(3), -1.)), .65);
    auto const identity = DenseMatrix::identity(3);
    auto const expected = add(add(scale(add(kron(adjoint(a), a), kron(a, adjoint(a))), -.7), kron(onsite, identity)),
                              kron(identity, onsite));
    for (auto scheme : schemes) {
        auto const encoding = mdl::make_qubit_encoding(model, scheme);
        EXPECT_LT(action_distance(pauli_matrix(mdl::export_pauli(model, encoding)), physical_indices(encoding), expected), 1e-13);
    }
}

TEST(PauliExport, TraceCompletionGroundEnergyAndDegenerateInvalidStates) {
    auto const site = mdl::spin(1.);
    mdl::Hamiltonian model(lat::chain(2), site);
    // Site 1 is unused: its invalid level can be tied at the ground energy.
    model.add(2. * L(site, "Sz", 0) * L(site, "Sz", 0) + .13 * L(site, "Sz", 0) + 4. * L(site, "I", 0));
    for (auto scheme : schemes) {
        auto const encoding = mdl::make_qubit_encoding(model, scheme);
        auto const matrix = pauli_matrix(mdl::export_pauli(model, encoding, exact_options()));
        auto const indices = physical_indices(encoding);
        std::vector<std::uint64_t> invalid;
        for (std::uint64_t j = 0; j < 16; ++j)
            if (std::find(indices.begin(), indices.end(), j) == indices.end()) invalid.push_back(j);
        EXPECT_NEAR(hermitian_eigenvalues(matrix).front(), 4., 1e-12);
        EXPECT_NEAR(hermitian_eigenvalues(restrict_to(matrix, invalid)).front(), 4., 1e-12);
        auto const penalty = pauli_matrix(mdl::invalid_state_penalty(encoding, .125));
        EXPECT_NEAR(hermitian_eigenvalues(restrict_to(add(matrix, penalty), invalid)).front(), 4.125, 1e-12);
        EXPECT_LT(max_column_norm(penalty, indices), 1e-14);
        // E(Sz^2) differs from E(Sz)^2 outside the physical space: products
        // belong in the physical algebra, before the local extension.
        auto const z = pauli_matrix(observable(model, L(site, "Sz", 0), encoding));
        auto const square = pauli_matrix(observable(model, L(site, "Sz", 0) * L(site, "Sz", 0), encoding));
        auto const mismatch = add(square, multiply(z, z), -1.);
        EXPECT_GT(frobenius_distance(mismatch, DenseMatrix(16)), .6);
        EXPECT_LT(max_column_norm(mismatch, indices), 1e-14);
    }
}

TEST(PauliExport, SpinTwoExactPruningWithTightWorkspace) {
    // Zero pruning and finite validation accuracy are different controls; a
    // tight workspace budget exercises uncached charge decompositions.
    auto const two = mdl::spin(2.);
    mdl::Hamiltonian model(lat::chain(2), two, mdl::FermionConvention::strict, {"Sz"});
    model.add(L(two, "Sz", 0) * L(two, "Sz", 1) + .5 * L(two, "Sp", 0) * L(two, "Sm", 1)
              + .5 * L(two, "Sm", 0) * L(two, "Sp", 1));
    auto options = exact_options();
    options.max_local_matrix_bytes = 8 * 8 * sizeof(Complex);
    auto const exported = mdl::export_pauli(model, options);
    EXPECT_EQ(exported.num_qubits, 6);
    EXPECT_NEAR(hermitian_eigenvalues(pauli_matrix(exported)).front(), -6., 1e-12);
}

TEST(PauliExport, AmplifiedNearIdentityFactorSurvives) {
    // Small local departures from identity must survive a large coefficient:
    // the transform prunes exact zeros only.
    auto const basis = mdl::physical_basis({mdl::spin()});
    mdl::LocalFactor nearly_identity{0, {1., 1e-14, 1e-14, 1.}};
    mdl::CompiledOperator const op{{2}, {{1e14, {nearly_identity}, {}}}};
    auto const precise = mdl::export_pauli(op, basis);
    ASSERT_EQ(precise.terms.count("X"), 1u);
    EXPECT_NEAR(std::abs(precise.terms.at("X") - Complex(1.)), 0., 1e-14);
    EXPECT_NEAR(precise.terms.at("I").real(), 1e14, 1.);
}

TEST(PauliExport, RepeatedSitesMultiplyInOrderAndDimensionsAreChecked) {
    auto const basis = mdl::physical_basis({mdl::spin(1.)});
    auto const sp = mdl::spin(1.)->op("Sp").matrix, sm = mdl::spin(1.)->op("Sm").matrix;
    mdl::CompiledOperator const product{{3}, {{1., {{0, sp}, {0, sm}}, {}}}};
    auto const [x, y, z] = spin_matrices(1.);
    auto const raising = add(x, y, Complex(0., 1.)), lowering = adjoint(raising);
    auto const expected = multiply(raising, lowering);
    QubitEncoding const encoding(basis);
    EXPECT_LT(action_distance(pauli_matrix(mdl::export_pauli(product, basis, encoding)), physical_indices(encoding), expected), 1e-13);
    expect_invalid([&] { (void)mdl::export_pauli(mdl::CompiledOperator{{3, 3}, {}}, basis); }, "site count");
    expect_invalid([&] { (void)mdl::export_pauli(mdl::CompiledOperator{{2}, {}}, basis); }, "local dimension");
    expect_invalid([&] { (void)mdl::export_pauli(mdl::CompiledOperator{{3}, {{1., {{1, sp}}, {}}}}, basis); }, "out of range");
    expect_invalid([&] { (void)mdl::export_pauli(mdl::CompiledOperator{{3}, {{1., {{0, {1., 0.}}}, {}}}}, basis); },
                   "dimension mismatch");
}
