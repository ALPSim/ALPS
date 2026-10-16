// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Classic molecular DMRG against exact references. The occupation-bit exact
// diagonalization and the RDM, trace, contraction and variational checks are
// adapted from skilledwolf/alpsx (private repository), branch
// WIP/classic-dmrg-chemistry at d4384a7f8, tests/dmrg/classic_molecular.cpp
// (a command-line program there, one ctest per case; GoogleTest cases here).
// Added here: the same Hamiltonian built from ALPS::models four-operator
// terms and solved by ALPS::ed, the validation cases, and the callback and
// scratch checks.
#include "support.hpp"
#include <gtest/gtest.h>

#include <alps/dmrg/molecular.hpp>
#include <alps/models/site.hpp>
#include <alps/models/terms.hpp>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace classic_test;

namespace {

// Set bits of an occupation word (C++17 has no std::popcount).
int popcount(unsigned bits) { return static_cast<int>(std::bitset<32>(bits).count()); }

struct Op { int mode; bool create; };
struct Monomial { double coefficient; std::vector<Op> factors; };

// Deterministic integrals with the real-orbital symmetries (a sum of
// separable pair channels, as in alpsx).
cd::MolecularInput input_for(int n, int na, int nb) {
    cd::MolecularInput input;
    input.orbitals = static_cast<std::size_t>(n);
    input.alpha_electrons = static_cast<std::size_t>(na);
    input.beta_electrons = static_cast<std::size_t>(nb);
    input.core_energy = .731;
    input.one_body.resize(static_cast<std::size_t>(n * n));
    input.two_body.resize(static_cast<std::size_t>(n * n * n * n));
    auto pair = [](int p, int q, int channel) { return .12 * std::cos((p + q + 1) * (channel + 1)) + .06 * (p == q); };
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            input.one_body[static_cast<std::size_t>(p * n + q)] = p == q ? .17 * (p + 1) * (p + 2) : -.41 / (1 + std::abs(p - q));
            for (int r = 0; r < n; ++r)
                for (int s = 0; s < n; ++s)
                    for (int channel = 0; channel < 3; ++channel)
                        input.two_body[static_cast<std::size_t>(((p * n + q) * n + r) * n + s)] +=
                            pair(p, q, channel) * pair(r, s, channel);
        }
    return input;
}

cd::Options exact(int kept) {
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(5, static_cast<std::size_t>(kept));
    options.warmup_states = static_cast<std::size_t>(kept);
    options.lanczos_tolerance = 1.e-12;
    return options;
}

// Occupation-bit ED on modes 2p (beta) and 2p+1 (alpha): the ground state and
// expectation values of normal-ordered monomials.
struct BitReference {
    std::vector<unsigned> states;
    std::map<unsigned, int> indices;
    Eigen::VectorXd psi;
    double energy = 0.;

    std::pair<int, double> apply(Monomial const& term, int col) const {
        unsigned bits = states[static_cast<std::size_t>(col)];
        double value = term.coefficient;
        for (auto op = term.factors.rbegin(); op != term.factors.rend(); ++op) {
            bool const occupied = (bits >> op->mode) & 1u;
            if (occupied == op->create) return {0, 0.};
            if (popcount(bits & ((1u << op->mode) - 1u)) % 2) value = -value;
            bits ^= 1u << op->mode;
        }
        return {indices.at(bits), value};
    }

    double expectation(Monomial const& term) const {
        double value = 0.;
        for (int col = 0; col < static_cast<int>(states.size()); ++col) {
            auto const [row, coefficient] = apply(term, col);
            value += psi[row] * coefficient * psi[col];
        }
        return value;
    }

    BitReference(cd::MolecularInput const& input) {
        int const n = static_cast<int>(input.orbitals);
        int const na = static_cast<int>(input.alpha_electrons), nb = static_cast<int>(input.beta_electrons);
        for (unsigned bits = 0; bits < (1u << (2 * n)); ++bits) {
            int a = 0, b = 0;
            for (int p = 0; p < n; ++p) {
                a += (bits >> (2 * p + 1)) & 1u;
                b += (bits >> (2 * p)) & 1u;
            }
            if (a == na && b == nb) {
                indices[bits] = static_cast<int>(states.size());
                states.push_back(bits);
            }
        }
        auto const size = static_cast<Eigen::Index>(states.size());
        Eigen::MatrixXd h = Eigen::MatrixXd::Zero(size, size);
        auto add = [&](Monomial const& term) {
            for (int col = 0; col < static_cast<int>(states.size()); ++col) {
                auto const [row, value] = apply(term, col);
                h(row, col) += value;
            }
        };
        auto at = [&](std::vector<double> const& v, int index) { return v[static_cast<std::size_t>(index)]; };
        for (int p = 0; p < n; ++p)
            for (int q = 0; q < n; ++q) {
                for (int a = 0; a < 2; ++a) add({at(input.one_body, p * n + q), {{2 * p + a, true}, {2 * q + a, false}}});
                for (int r = 0; r < n; ++r)
                    for (int s = 0; s < n; ++s)
                        for (int a = 0; a < 2; ++a)
                            for (int b = 0; b < 2; ++b)
                                add({.5 * at(input.two_body, ((p * n + q) * n + r) * n + s),
                                     {{2 * p + a, true}, {2 * r + b, true}, {2 * s + b, false}, {2 * q + a, false}}});
            }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ed(h);
        if (ed.info() != Eigen::Success) throw std::runtime_error("ED reference did not converge");
        psi = ed.eigenvectors().col(0);
        energy = ed.eigenvalues()[0] + input.core_energy;
    }
};

struct Errors {
    double energy, rdm, rdm_energy, trace_alpha, trace_beta, contraction;
};

Errors compare(cd::MolecularInput const& input, cd::MolecularResult const& result, BitReference const& reference) {
    int const n = static_cast<int>(input.orbitals);
    int const ne = static_cast<int>(input.alpha_electrons + input.beta_electrons);
    auto at = [](std::vector<double> const& v, int index) { return v[static_cast<std::size_t>(index)]; };
    Errors errors{std::abs(result.energy - reference.energy), 0., 0., 0., 0., 0.};
    double from_rdms = input.core_energy;
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            double const aa = reference.expectation({1., {{2 * q + 1, true}, {2 * p + 1, false}}});
            double const bb = reference.expectation({1., {{2 * q, true}, {2 * p, false}}});
            errors.rdm = std::max({errors.rdm, std::abs(aa - at(result.one_rdm_alpha, p * n + q)),
                                   std::abs(bb - at(result.one_rdm_beta, p * n + q))});
            from_rdms += at(input.one_body, p * n + q) * (at(result.one_rdm_alpha, q * n + p) + at(result.one_rdm_beta, q * n + p));
            for (int r = 0; r < n; ++r)
                for (int s = 0; s < n; ++s) {
                    double expected = 0.;
                    for (int a = 0; a < 2; ++a)
                        for (int b = 0; b < 2; ++b)
                            expected += reference.expectation(
                                {1., {{2 * p + a, true}, {2 * r + b, true}, {2 * s + b, false}, {2 * q + a, false}}});
                    int const index = ((p * n + q) * n + r) * n + s;
                    errors.rdm = std::max(errors.rdm, std::abs(expected - at(result.two_rdm, index)));
                    from_rdms += .5 * at(input.two_body, index) * at(result.two_rdm, index);
                }
        }
    errors.rdm_energy = std::abs(result.energy - from_rdms);
    double trace_a = 0., trace_b = 0.;
    for (int p = 0; p < n; ++p) {
        trace_a += at(result.one_rdm_alpha, p * n + p);
        trace_b += at(result.one_rdm_beta, p * n + p);
        for (int q = 0; q < n; ++q) {
            double contracted = 0.;
            for (int r = 0; r < n; ++r) contracted += at(result.two_rdm, ((p * n + q) * n + r) * n + r);
            double const one = at(result.one_rdm_alpha, q * n + p) + at(result.one_rdm_beta, q * n + p);
            errors.contraction = std::max(errors.contraction, std::abs(contracted - (ne - 1) * one));
        }
    }
    errors.trace_alpha = std::abs(trace_a - static_cast<double>(input.alpha_electrons));
    errors.trace_beta = std::abs(trace_b - static_cast<double>(input.beta_electrons));
    return errors;
}

void expect_exact(int n, int na, int nb) {
    auto const input = input_for(n, na, nb);
    BitReference const reference(input);
    auto const result = cd::solve_molecular(input, exact(128));
    ASSERT_EQ(result.orbitals, static_cast<std::size_t>(n));
    ASSERT_EQ(result.two_rdm.size(), static_cast<std::size_t>(n * n * n * n));
    auto const errors = compare(input, result, reference);
    EXPECT_LT(errors.energy, 1e-9) << n << " orbitals (" << na << "," << nb << ")";
    EXPECT_LT(errors.rdm, 1e-8);
    EXPECT_LT(errors.rdm_energy, 1e-9);
    EXPECT_LT(errors.trace_alpha, 1e-9);
    EXPECT_LT(errors.trace_beta, 1e-9);
    EXPECT_LT(errors.contraction, 1e-8);
    EXPECT_FALSE(result.iterations.empty());
}

// Added here: the spin blocks of the 2-RDM, the caller's observables and the
// CI vector against the same exact reference. Modes are 2p (down) and 2p+1
// (up), the ALPS order the CI vector is documented in.
void expect_parity(int n, int na, int nb) {
    auto input = input_for(n, na, nb);
    auto const fs = mdl::spinful_fermion();
    auto L = [&](char const* name, int site) { return mdl::TermSum::located(fs, name, site); };
    int const last = n - 1;
    input.observables["szsz"] = L("Sz", 0) * L("Sz", last);
    input.observables["doublons"] = L("n_updn", 0) * L("n_updn", last) * 2.;
    input.observables["hopping"] = L("c_dag_up", 0) * L("c_up", last) + L("c_dag_up", last) * L("c_up", 0);
    input.observables["charge_changing"] = L("c_dag_up", 0);
    if (n >= 4) input.observables["four_sites"] = L("n_up", 0) * L("n_dn", 1) * L("Sp", 2) * L("Sm", 3);
    input.ci_vector = true;
    BitReference const reference(input);
    auto const result = cd::solve_molecular(input, exact(128));
    auto at = [](std::vector<double> const& v, int index) { return v[static_cast<std::size_t>(index)]; };
    double blocks = 0.;
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q)
            for (int r = 0; r < n; ++r)
                for (int s = 0; s < n; ++s) {
                    int const index = ((p * n + q) * n + r) * n + s;
                    auto element = [&](int a, int b) {
                        return reference.expectation({1., {{2 * p + a, true}, {2 * r + b, true}, {2 * s + b, false}, {2 * q + a, false}}});
                    };
                    blocks = std::max({blocks, std::abs(element(1, 1) - at(result.two_rdm_aa, index)),
                                       std::abs(element(1, 0) - at(result.two_rdm_ab, index)),
                                       std::abs(element(0, 0) - at(result.two_rdm_bb, index)),
                                       std::abs(at(result.two_rdm_aa, index) + at(result.two_rdm_ab, index) +
                                                at(result.two_rdm_ab, ((r * n + s) * n + p) * n + q) +
                                                at(result.two_rdm_bb, index) - at(result.two_rdm, index))});
                }
    EXPECT_LT(blocks, 1e-8) << n << " orbitals (" << na << "," << nb << ")";
    auto number = [](int mode) { return std::vector<Op>{{mode, true}, {mode, false}}; };
    auto product = [](std::vector<std::vector<Op>> const& parts, double coefficient) {
        Monomial term{coefficient, {}};
        for (auto const& part : parts) term.factors.insert(term.factors.end(), part.begin(), part.end());
        return term;
    };
    double szsz = 0.;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            szsz += reference.expectation(product({number(2 * 0 + a), number(2 * last + b)}, .25 * (a ? 1 : -1) * (b ? 1 : -1)));
    EXPECT_NEAR(result.observables.at("szsz"), szsz, 1e-8);
    EXPECT_NEAR(result.observables.at("doublons"),
                reference.expectation(product({number(1), number(0), number(2 * last + 1), number(2 * last)}, 2.)), 1e-8);
    EXPECT_NEAR(result.observables.at("hopping"),
                reference.expectation({1., {{1, true}, {2 * last + 1, false}}}) +
                    reference.expectation({1., {{2 * last + 1, true}, {1, false}}}), 1e-8);
    EXPECT_EQ(result.observables.at("charge_changing"), 0.);
    if (n >= 4)
        EXPECT_NEAR(result.observables.at("four_sites"),
                    reference.expectation(product({number(1), number(2), {{5, true}, {4, false}}, {{6, true}, {7, false}}}, 1.)),
                    1e-8);
    // The CI vector, sign by sign, up to the reference's global sign.
    ASSERT_EQ(result.ci_coefficients.size(), reference.states.size());
    int largest = 0;
    for (int i = 0; i < static_cast<int>(reference.psi.size()); ++i)
        if (std::abs(reference.psi[i]) > std::abs(reference.psi[largest])) largest = i;
    double const sign = reference.psi[largest] < 0. ? -1. : 1.;
    double ci = 0.;
    for (std::size_t k = 0; k < result.ci_coefficients.size(); ++k) {
        unsigned bits = 0;
        for (int p = 0; p < n; ++p)
            bits |= static_cast<unsigned>(((result.ci_alpha[k] >> p) & 1u) << (2 * p + 1)) |
                    static_cast<unsigned>(((result.ci_beta[k] >> p) & 1u) << (2 * p));
        if (k) EXPECT_TRUE(result.ci_alpha[k - 1] < result.ci_alpha[k] ||
                           (result.ci_alpha[k - 1] == result.ci_alpha[k] && result.ci_beta[k - 1] < result.ci_beta[k]));
        ci = std::max(ci, std::abs(result.ci_coefficients[k] - sign * reference.psi[reference.indices.at(bits)]));
    }
    EXPECT_LT(ci, 1e-8) << n << " orbitals (" << na << "," << nb << ")";
}

} // namespace

TEST(ClassicMolecular, SpinBlocksObservablesAndCiVectorMatchExactDiagonalization) {
    expect_parity(1, 1, 0);
    expect_parity(2, 1, 1);
    expect_parity(3, 2, 1);
    expect_parity(4, 2, 2);
    expect_parity(5, 3, 2);
    expect_parity(6, 3, 3);
}

TEST(ClassicMolecular, TruncatedCiVectorCarriesTheRunsOwnRdms) {
    // With 8 kept states the state differs from the exact one; the CI vector
    // must still be the state the RDMs were measured on.
    auto input = input_for(6, 3, 2);
    input.ci_vector = true;
    BitReference reference(input);
    double const exact_energy = reference.energy;
    auto const result = cd::solve_molecular(input, exact(8));
    EXPECT_GT(result.energy - exact_energy, 1e-6);
    ASSERT_EQ(result.ci_coefficients.size(), reference.states.size());
    double norm = 0.;
    for (std::size_t k = 0; k < result.ci_coefficients.size(); ++k) {
        unsigned bits = 0;
        for (int p = 0; p < 6; ++p)
            bits |= static_cast<unsigned>(((result.ci_alpha[k] >> p) & 1u) << (2 * p + 1)) |
                    static_cast<unsigned>(((result.ci_beta[k] >> p) & 1u) << (2 * p));
        reference.psi[reference.indices.at(bits)] = result.ci_coefficients[k];
        norm += result.ci_coefficients[k] * result.ci_coefficients[k];
    }
    EXPECT_NEAR(norm, 1., 1e-12);
    double error = 0.;
    for (int p = 0; p < 6; ++p)
        for (int q = 0; q < 6; ++q) {
            error = std::max(error, std::abs(reference.expectation({1., {{2 * q + 1, true}, {2 * p + 1, false}}}) -
                                             result.one_rdm_alpha[static_cast<std::size_t>(p * 6 + q)]));
            error = std::max(error, std::abs(reference.expectation({1., {{2 * q, true}, {2 * p, false}}}) -
                                             result.one_rdm_beta[static_cast<std::size_t>(p * 6 + q)]));
        }
    EXPECT_LT(error, 1e-8);
}

TEST(ClassicMolecular, RejectedObservables) {
    auto const fs = mdl::spinful_fermion();
    auto L = [&](char const* name, int site) { return mdl::TermSum::located(fs, name, site); };
    auto rejects = [](mdl::TermSum const& observable) {
        auto input = input_for(5, 2, 2);
        input.observables["bad"] = observable;
        EXPECT_THROW(cd::solve_molecular(input, exact(16)), std::invalid_argument);
    };
    rejects(L("n", 5));                                                      // outside the orbitals
    rejects(L("n", 0) * L("n", 1) * L("n", 2) * L("n", 3) * L("n", 4));      // five sites
    rejects(L("n", 0) * std::complex<double>(0., 1.));                      // complex coefficient
    rejects(mdl::TermSum::located(mdl::spin(0.5), "Sz", 0));                // foreign local basis
}

TEST(ClassicMolecular, OneOrbitalWithPadding) {
    expect_exact(1, 1, 1);
    expect_exact(1, 0, 0);
}

TEST(ClassicMolecular, SmallSpacesMatchExactDiagonalization) {
    expect_exact(2, 1, 1);
    expect_exact(3, 2, 1);
    expect_exact(4, 2, 1);
}

TEST(ClassicMolecular, FiveOrbitalsMatchExactDiagonalization) { expect_exact(5, 2, 2); }

TEST(ClassicMolecular, TruncatedRunStaysVariationalAndConsistent) {
    auto const input = input_for(5, 2, 2);
    BitReference const reference(input);
    auto const result = cd::solve_molecular(input, exact(8));
    auto const errors = compare(input, result, reference);
    EXPECT_GE(result.energy - reference.energy, -1e-9);
    EXPECT_LT(errors.rdm_energy, 1e-9);
    EXPECT_LT(errors.trace_alpha, 1e-9);
    EXPECT_LT(errors.trace_beta, 1e-9);
    EXPECT_LT(errors.contraction, 1e-8);
}

TEST(ClassicMolecular, AgreesWithTheModelsHamiltonianSolvedByEd) {
    // The same operator from ALPS::models four-operator terms and ALPS::ed.
    int const n = 4;
    auto const input = input_for(n, 2, 2);
    auto const fs = mdl::spinful_fermion();
    mdl::Hamiltonian model(lat::chain(n, 1.0, lat::Boundary::open), fs, mdl::FermionConvention::strict, {"N", "Sz"});
    char const* create[2] = {"c_dag_up", "c_dag_dn"};
    char const* destroy[2] = {"c_up", "c_dn"};
    auto L = [&](char const* name, int site) { return mdl::TermSum::located(fs, name, site); };
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            for (int a = 0; a < 2; ++a)
                model.add(L(create[a], p) * L(destroy[a], q) * input.one_body[static_cast<std::size_t>(p * n + q)]);
            for (int r = 0; r < n; ++r)
                for (int s = 0; s < n; ++s) {
                    double const v = input.two_body[static_cast<std::size_t>(((p * n + q) * n + r) * n + s)];
                    for (int a = 0; a < 2; ++a)
                        for (int b = 0; b < 2; ++b)
                            model.add(L(create[a], p) * L(create[b], r) * L(destroy[b], s) * L(destroy[a], q) * (.5 * v));
                }
        }
    auto const ed = ed_energies(model, {{"N", 4.}, {"Sz", 0.}});
    auto const result = cd::solve_molecular(input, exact(128));
    EXPECT_NEAR(result.energy - input.core_energy, ed.at(0), 1e-9);
}

TEST(ClassicMolecular, CallbackObservesAndAbortsAndScratchIsRemoved) {
    auto const input = input_for(4, 2, 1);
    auto const parent = std::filesystem::temp_directory_path() / "alps-classic-molecular-test";
    std::filesystem::remove_all(parent);
    std::filesystem::create_directories(parent);
    auto options = exact(32);
    options.scratch_parent = parent;
    std::size_t seen = 0;
    auto const result = cd::solve_molecular(input, options, [&](cd::Iteration const&) { ++seen; });
    EXPECT_EQ(seen, result.iterations.size());
    EXPECT_TRUE(std::filesystem::is_empty(parent));
    struct Stop {};
    EXPECT_THROW(cd::solve_molecular(input, options, [](cd::Iteration const& it) { if (it.index >= 2) throw Stop{}; }), Stop);
    EXPECT_TRUE(std::filesystem::is_empty(parent));
    // The engine state is restored: the next solve runs normally.
    EXPECT_NEAR(cd::solve_molecular(input, exact(128)).energy, BitReference(input).energy, 1e-9);
    std::filesystem::remove_all(parent);
}

TEST(ClassicMolecular, RejectedInputs) {
    auto const good = input_for(3, 1, 1);
    auto rejects = [](cd::MolecularInput const& input, cd::Options const& options = exact(16)) {
        EXPECT_THROW(cd::solve_molecular(input, options), std::invalid_argument);
    };
    auto input = good;
    input.one_body[1] += 1e-3; // h[0,1] != h[1,0]
    rejects(input);
    input = good;
    input.two_body[1] += 1e-3; // (00|01) != (00|10)
    rejects(input);
    input = good;
    input.alpha_electrons = 4;
    rejects(input);
    input = good;
    input.two_body.pop_back();
    rejects(input);
    input = good;
    input.one_body[0] = std::numeric_limits<double>::quiet_NaN();
    rejects(input);
    input = good;
    input.orbitals = 0;
    input.one_body.clear();
    input.two_body.clear();
    rejects(input);
    auto options = exact(16);
    options.sector = {{"N", 2.}};
    rejects(good, options);
    options = exact(16);
    options.num_eigenvalues = 2;
    rejects(good, options);
    options = exact(16);
    options.arithmetic = cd::Arithmetic::complex;
    rejects(good, options);
    options = exact(16);
    options.truncation_max_states = 8;
    rejects(good, options);
    // An asymmetric roundoff entry is rejected unless the explicit cutoff
    // prunes it.
    input = good;
    std::fill(input.two_body.begin(), input.two_body.end(), 0.);
    input.two_body[1] = 1e-14;
    rejects(input);
    input.integral_cutoff = 1e-12;
    EXPECT_NO_THROW(cd::solve_molecular(input, exact(16)));
}
