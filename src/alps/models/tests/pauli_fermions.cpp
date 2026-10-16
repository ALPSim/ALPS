// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Occupation-bit references for the fermionic Pauli export: complex hopping,
// pairing, a three-body term and local phases of spinless and spinful chains
// against the canonical anticommutation relations applied directly to
// occupied states (no operator compiler, no fermion-to-qubit library), odd
// observables, the register layout, validation and the absence of a global
// intermediate. Cases follow skilledwolf/alpsx (branch
// agent/v3-execution-and-api-audit at f8ee6e98b,
// src/model/tests/models/pauli_export.cpp and tests/alps/test_fermion_pauli.py),
// adapted to GoogleTest and the CompiledOperator/PhysicalBasis entry.
#include "pauli_support.hpp"
#include <alps/models/builtin.hpp>
#include <alps/models/hamiltonian.hpp>
#include <alps/models/pauli.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;
using namespace alps_test;
using mdl::QubitEncoding;
using mdl::QubitEncodingScheme;

namespace {
mdl::TermSum L(mdl::SiteBasisPtr const& site, char const* name, int i) { return mdl::TermSum::located(site, name, i); }

mdl::PauliExportOptions observe() {
    mdl::PauliExportOptions options;
    options.require_hermitian = options.require_even = options.require_conservation = false;
    options.tolerance = 0.;
    return options;
}
mdl::PauliOperator observable(mdl::Hamiltonian const& model, mdl::TermSum const& op) {
    return mdl::export_pauli(model.operator_decomposition(op).op, mdl::physical_basis(model), observe());
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

struct Reference {
    mdl::Hamiltonian model;
    DenseMatrix matrix;
    std::vector<DenseMatrix> annihilators; // CAR order
    std::vector<mdl::TermSum> operators;   // the kit's annihilators, same order
};

// Open four-site spinless chain: site s is qubit s.
Reference spinless() {
    auto const site = mdl::spinless_fermion();
    mdl::Hamiltonian model(lat::chain(4, 1.0, lat::Boundary::open), site);
    std::vector<DenseMatrix> c;
    for (int i = 0; i < 4; ++i) c.push_back(annihilator(4, i));
    DenseMatrix reference(16);
    struct Hop { int i, j; Complex t; };
    for (auto const& [i, j, t] : {Hop{0, 1, .7}, Hop{1, 3, Complex(-.3, .21)}, Hop{3, 0, Complex(0., .11)}}) {
        model.add(t * L(site, "c_dag", i) * L(site, "c", j), true);
        auto const term = scale(multiply(adjoint(c[static_cast<std::size_t>(i)]), c[static_cast<std::size_t>(j)]), t);
        reference = add(add(reference, term), adjoint(term));
    }
    model.add(Complex(.23, .17) * L(site, "c", 3) * L(site, "c", 0), true);
    auto const pairing = scale(multiply(c[3], c[0]), Complex(.23, .17));
    reference = add(add(reference, pairing), adjoint(pairing));
    auto const number = [&](std::size_t i) { return multiply(adjoint(c[i]), c[i]); };
    model.add(.41 * L(site, "n", 0) * L(site, "n", 2) * L(site, "n", 3));
    reference = add(reference, multiply(multiply(number(0), number(2)), number(3)), .41);
    model.add(.13 * L(site, "I", 0));
    reference = add(reference, DenseMatrix::identity(16), .13);
    std::vector<mdl::TermSum> operators;
    for (int i = 0; i < 4; ++i) operators.push_back(L(site, "c", i));
    return {std::move(model), std::move(reference), std::move(c), std::move(operators)};
}

// Open two-site spinful chain. Local digits 0,1,2,3 = empty, up, down,
// double: the up bit is the lower qubit of each site; site-major CAR order
// is (down, up), i.e. modes (dn0, up0, dn1, up1) on qubits (1, 0, 3, 2).
Reference spinful() {
    auto const site = mdl::spinful_fermion();
    mdl::Hamiltonian model(lat::chain(2, 1.0, lat::Boundary::open), site);
    std::vector<DenseMatrix> c;
    for (int m = 0; m < 4; ++m) c.push_back(annihilator(4, m, {1, 0, 3, 2}));
    model.add(Complex(.37, .19) * L(site, "c_dag_up", 0) * L(site, "c_dn", 1), true);
    model.add(Complex(0., .29) * L(site, "c_dag_dn", 0) * L(site, "c_up", 0), true);
    model.add(.17 * L(site, "c_dn", 1) * L(site, "c_up", 0), true);
    auto term = scale(multiply(adjoint(c[1]), c[2]), Complex(.37, .19));
    term = add(term, scale(multiply(adjoint(c[0]), c[1]), Complex(0., .29)));
    term = add(term, scale(multiply(c[2], c[1]), .17));
    auto const reference = add(term, adjoint(term));
    std::vector<mdl::TermSum> operators{L(site, "c_dn", 0), L(site, "c_up", 0), L(site, "c_dn", 1), L(site, "c_up", 1)};
    return {std::move(model), reference, std::move(c), std::move(operators)};
}

void check_car(std::vector<DenseMatrix> const& a, double tolerance) {
    auto const identity = DenseMatrix::identity(a.front().dim);
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a.size(); ++j) {
            EXPECT_LT(max_abs(add(multiply(a[i], a[j]), multiply(a[j], a[i]))), tolerance) << i << " " << j;
            auto const anticommutator = add(multiply(a[i], adjoint(a[j])), multiply(adjoint(a[j]), a[i]));
            EXPECT_LT(max_abs_difference(anticommutator, i == j ? identity : DenseMatrix(identity.dim)), tolerance) << i << " " << j;
        }
}
} // namespace

TEST(PauliFermions, ReferencesSatisfyTheCanonicalAnticommutationRelations) {
    check_car(spinless().annihilators, 1e-14);
    check_car(spinful().annihilators, 1e-14);
    // A spectrum would not reveal a basis error; the full action does.
    auto const reference = spinless().matrix;
    DenseMatrix reversed(16);
    for (std::size_t r = 0; r < 16; ++r)
        for (std::size_t c = 0; c < 16; ++c) reversed(r, c) = reference(15 - r, 15 - c);
    EXPECT_GT(max_abs_difference(reference, reversed), .1);
}

TEST(PauliFermions, SpinlessAndSpinfulChainsActLikeTheOccupationReferences) {
    for (auto const& make : {&spinless, &spinful}) {
        auto const reference = make();
        auto const exported = mdl::export_pauli(reference.model);
        EXPECT_EQ(exported.num_qubits, 4);
        // Every code word of a dim-2/4 occupation basis is physical.
        EXPECT_LT(max_abs_difference(pauli_matrix(exported), reference.matrix), 2e-14);
        EXPECT_LT(exported.discarded_l1, 1e-13);
        std::vector<DenseMatrix> annihilators;
        for (auto const& op : reference.operators) annihilators.push_back(pauli_matrix(observable(reference.model, op)));
        for (std::size_t m = 0; m < annihilators.size(); ++m)
            EXPECT_LT(max_abs_difference(annihilators[m], reference.annihilators[m]), 2e-14) << "mode " << m;
        check_car(annihilators, 2e-14);
    }
    EXPECT_EQ(mdl::export_pauli(spinless().model).state_index({1, 0, 1, 0}), 5u);
    EXPECT_EQ(mdl::export_pauli(spinful().model).fermion_mode_qubits, (std::vector<int>{1, 0, 3, 2}));
}

TEST(PauliFermions, HubbardSpectrumMatchesTheDenseAssembly) {
    auto const model = mdl::hubbard(lat::chain(3, 1.0, lat::Boundary::open), 1., 4., .7);
    auto const compiled = mdl::compile_hamiltonian(model);
    auto const reference = hermitian_eigenvalues(dense_hamiltonian(model, compiled));
    auto const exported = mdl::export_pauli(model);
    EXPECT_EQ(exported.num_qubits, 6);
    auto const encoding = mdl::make_qubit_encoding(model);
    auto const spectrum = hermitian_eigenvalues(restrict_to(pauli_matrix(exported), physical_indices(encoding)));
    ASSERT_EQ(spectrum.size(), reference.size());
    for (std::size_t k = 0; k < spectrum.size(); ++k) EXPECT_NEAR(spectrum[k], reference[k], 1e-12);
}

TEST(PauliFermions, ValidationAndNoDenseIntermediate) {
    auto const site = mdl::spinless_fermion();
    auto const reference = spinless();
    auto const compiled = mdl::compile_hamiltonian(reference.model);
    // Pairing breaks a stale particle-number declaration.
    auto stale = mdl::physical_basis(reference.model);
    stale.conserved_quantum_numbers = {"N"};
    expect_invalid([&] { (void)mdl::export_pauli(compiled.op, stale); }, "does not conserve declared quantum number N");
    mdl::Hamiltonian declared(lat::chain(4, 1.0, lat::Boundary::open), site, mdl::FermionConvention::strict, {"N"});
    declared.add(reference.model.expand());
    EXPECT_THROW(mdl::export_pauli(declared), std::invalid_argument);
    auto lattice_convention = mdl::physical_basis(reference.model);
    lattice_convention.fermion_convention = mdl::FermionConvention::lattice;
    expect_invalid([&] { (void)mdl::export_pauli(compiled.op, lattice_convention); }, "strict");
    std::vector<mdl::SiteBasisPtr> const pair{site, site};
    auto const basis = mdl::physical_basis(pair);
    auto const odd = mdl::compile_operator(L(site, "c", 0) + L(site, "c_dag", 0), pair);
    expect_invalid([&] { (void)mdl::export_pauli(odd, basis); }, "even fermionic");
    auto const skew = mdl::compile_operator(Complex(0., 1.) * L(site, "n", 0), pair);
    expect_invalid([&] { (void)mdl::export_pauli(skew, basis); }, "Hermitian");
    for (double bad : {-1., std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        mdl::PauliExportOptions options;
        options.tolerance = bad;
        expect_invalid([&] { (void)mdl::export_pauli(reference.model, options); }, "tolerance");
    }
    expect_invalid([&] { (void)mdl::make_qubit_encoding(reference.model, QubitEncodingScheme::gray); }, "binary dim-2/4");

    // A 128-site hopping chain yields 254 strings from three local decompositions.
    mdl::Hamiltonian large(lat::chain(128, 1.0, lat::Boundary::open), site);
    large.add_bond(-1., [site](int i, int j) { return L(site, "c_dag", i) * L(site, "c", j); }, std::nullopt, true);
    auto const exported = mdl::export_pauli(large);
    EXPECT_EQ(exported.num_qubits, 128);
    EXPECT_EQ(exported.terms.size(), 254u);
    EXPECT_LE(exported.local_decompositions, 3u);
    mdl::PauliExportOptions bounded;
    bounded.max_terms = 3;
    EXPECT_THROW(mdl::export_pauli(large, bounded), std::length_error);
}

TEST(PauliFermions, ZeroConstantsPruningAndLayout) {
    auto const site = mdl::spinful_fermion();
    mdl::Hamiltonian model(lat::chain(2), site);
    auto const empty = mdl::export_pauli(model);
    EXPECT_EQ(empty.num_qubits, 4);
    EXPECT_TRUE(empty.terms.empty());
    EXPECT_EQ(empty.site_qubits, (std::vector<int>{2, 2}));
    EXPECT_EQ(empty.fermion_mode_qubits, (std::vector<int>{1, 0, 3, 2}));
    EXPECT_EQ(empty.state_index({1, 2}), 9u);
    EXPECT_THROW(empty.state_index({4, 0}), std::invalid_argument);
    model.add(2. * L(site, "I", 0));
    model.add(1e-13 * L(site, "n", 0));
    auto const result = mdl::export_pauli(model);
    ASSERT_EQ(result.terms.size(), 1u);
    EXPECT_NEAR(result.terms.at("IIII").real(), 2., 1e-12);
    EXPECT_NEAR(result.discarded_l1, 1e-13, 1e-15);
    EXPECT_NEAR(result.coefficient_l1(), 2., 1e-12);
    EXPECT_EQ(result.coefficient_l1(false), 0.);
}

TEST(PauliFermions, RejectsNonFiniteAndMixedParity) {
    auto const site = mdl::spinless_fermion();
    std::vector<mdl::SiteBasisPtr> const pair{site, site};
    auto const basis = mdl::physical_basis(pair);
    auto const n = site->op("n").matrix;
    for (double value : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        mdl::CompiledOperator const op{{2, 2}, {{value, {{0, n}}, {}}}};
        expect_invalid([&] { (void)mdl::export_pauli(op, basis); }, "non-finite");
        mdl::CompiledOperator const entries{{2, 2}, {{1., {{0, {value, 0., 0., 1.}}}, {}}}};
        expect_invalid([&] { (void)mdl::export_pauli(entries, basis); }, "non-finite");
        mdl::Hamiltonian model(lat::chain(2), site);
        model.add(value * L(site, "n", 0));
        EXPECT_THROW(mdl::export_pauli(model), std::invalid_argument);
    }
    // An entry connecting both parities has no meaning under the CAR.
    mdl::CompiledOperator const mixed{{2, 2}, {{1., {{0, {1., 1., 0., 0.}}}, {}}}};
    expect_invalid([&] { (void)mdl::export_pauli(mixed, basis, observe()); }, "mixes fermionic parities");
    mdl::CompiledOperator const string{{2, 2}, {{1., {{0, n}}, {2}}}};
    expect_invalid([&] { (void)mdl::export_pauli(string, basis); }, "Jordan-Wigner site");
}
