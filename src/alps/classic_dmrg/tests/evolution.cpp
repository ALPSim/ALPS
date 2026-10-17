// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Time evolution of the classic DMRG against dense references on the ED
// basis of the initial state's sector. With enough kept states the sweeps
// truncate nothing, so every evolved state must be the second-order Trotter
// product exp(-dt/2 A) exp(-dt B) exp(-dt/2 A) applied n times (A the even
// bonds, B the odd ones, one-site terms split between a site's bonds as the
// gates split them) to rounding, and the exact exp(-i H t) to the Trotter
// error, which shrinks as dt^2. Spin chains (an XXZ quench, the XX domain
// wall), the Hubbard chain and spinless fermions (the gates' parity
// strings), drives, imaginary time, the exported states (overlaps and the
// Loschmidt amplitude), the schedule's bookkeeping and the rejections.
#include "support.hpp"
#include <gtest/gtest.h>

#include <alps/models/compiled.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <set>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

using namespace classic_test;

namespace {

using Complex = std::complex<double>;
using Dense = Eigen::MatrixXcd;
using Vector = Eigen::VectorXcd;

Dense dense(alps::ed::SparseMatrix<Complex> const& m) {
    Dense out = Dense::Zero(static_cast<Eigen::Index>(m.rows), static_cast<Eigen::Index>(m.columns));
    for (std::size_t r = 0; r < m.rows; ++r)
        for (auto k = m.offsets[r]; k < m.offsets[r + 1]; ++k)
            out(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(m.indices[static_cast<std::size_t>(k)])) +=
                m.values[static_cast<std::size_t>(k)];
    return out;
}

// exp(factor * h) of a Hermitian h.
Dense exponential(Dense const& h, Complex factor) {
    Eigen::SelfAdjointEigenSolver<Dense> solver(h);
    Vector phases(h.rows());
    for (Eigen::Index i = 0; i < h.rows(); ++i) phases[i] = std::exp(factor * solver.eigenvalues()[i]);
    return solver.eigenvectors() * phases.asDiagonal() * solver.eigenvectors().adjoint();
}

// The terms of `h` on the even and on the odd bonds, a one-site term divided
// evenly between the bonds of its site.
std::pair<mdl::CompiledOperator, mdl::CompiledOperator> even_odd(mdl::CompiledOperator const& h) {
    int const length = static_cast<int>(h.site_dims.size());
    mdl::CompiledOperator even{h.site_dims, {}}, odd{h.site_dims, {}};
    auto add = [&](int bond, mdl::LocalTerm term) { (bond % 2 == 0 ? even : odd).terms.push_back(std::move(term)); };
    for (auto const& term : h.terms) {
        if (term.factors.empty()) continue; // a phase
        if (term.factors.size() == 2) {
            add(std::min(term.factors[0].site, term.factors[1].site), term);
            continue;
        }
        int const site = term.factors[0].site;
        std::vector<int> bonds;
        if (site > 0) bonds.push_back(site - 1);
        if (site + 1 < length) bonds.push_back(site);
        for (int bond : bonds) {
            auto part = term;
            part.coefficient /= static_cast<double>(bonds.size());
            add(bond, part);
        }
    }
    return {even, odd};
}

mdl::CompiledOperator scaled(mdl::CompiledOperator op, double factor) {
    for (auto& term : op.terms) term.coefficient *= factor;
    return op;
}

mdl::CompiledOperator sum(mdl::CompiledOperator a, mdl::CompiledOperator const& b) {
    a.terms.insert(a.terms.end(), b.terms.begin(), b.terms.end());
    return a;
}

// The initial state's sector basis with dense references on it.
struct Reference {
    alps::ed::ProductEigenResult states;
    Vector initial;

    Reference(mdl::CompiledOperator const& h0, mdl::PhysicalBasis const& basis, Sector const& sector)
        : states(prepare(h0, basis, sector)) {
        std::size_t const dim = states.problem.dimension();
        initial.resize(static_cast<Eigen::Index>(dim));
        std::visit([&](auto const& values) {
            for (std::size_t i = 0; i < dim; ++i) initial[static_cast<Eigen::Index>(i)] = values[i];
        }, states.eigenvectors);
    }

    static alps::ed::ProductEigenResult prepare(mdl::CompiledOperator const& h0, mdl::PhysicalBasis const& basis,
                                                Sector const& sector) {
        alps::ed::ProductOptions options;
        options.sector = sector;
        alps::ed::EigenstateOptions request;
        request.method = alps::ed::EigenstateMethod::dense;
        return alps::ed::prepare_product(h0, basis, options).eigenstates(request);
    }

    Dense matrix(mdl::CompiledOperator const& op) const { return dense(states.problem.operator_matrix(op)); }

    double value(Vector const& psi, mdl::CompiledOperator const& op) const {
        return (psi.adjoint() * matrix(op) * psi)(0, 0).real() / psi.squaredNorm();
    }

    // One value per term, as Measurements::terms reports them.
    std::vector<double> values(Vector const& psi, mdl::CompiledOperator const& op) const {
        std::vector<double> out;
        for (auto const& term : op.terms) out.push_back(value(psi, mdl::CompiledOperator{op.site_dims, {term}}));
        return out;
    }

    // An exported state on this sector's basis.
    Vector vector(cd::MatrixProductState const& mps) const;
};

// The state on the whole product basis, site 0 most significant.
std::vector<Complex> contract(cd::MatrixProductState const& mps) {
    std::vector<Complex> carry{1.};
    std::size_t bond = 1;
    for (auto const& site : mps.sites) {
        std::vector<Complex> t = site.complex;
        if (t.empty()) t.assign(site.real.begin(), site.real.end());
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
    return carry;
}

Vector Reference::vector(cd::MatrixProductState const& mps) const {
    auto const full = contract(mps);
    auto const& dims = states.problem.site_dims();
    auto const& configurations = states.problem.states();
    std::size_t const length = dims.size(), dim = states.problem.dimension();
    Vector out(static_cast<Eigen::Index>(dim));
    double inside = 0., total = 0.;
    for (auto const& x : full) total += std::norm(x);
    for (std::size_t k = 0; k < dim; ++k) {
        std::size_t index = 0;
        for (std::size_t i = 0; i < length; ++i)
            index = index * static_cast<std::size_t>(dims[i]) + static_cast<std::size_t>(configurations[k * length + i]);
        out[static_cast<Eigen::Index>(k)] = full[index];
        inside += std::norm(full[index]);
    }
    EXPECT_NEAR(inside, total, 1e-12) << "the exported state leaves the sector";
    return out;
}

// A model's compiled Hamiltonian and its basis.
struct Model {
    mdl::CompiledOperator op;
    mdl::PhysicalBasis basis;
    explicit Model(mdl::Hamiltonian const& h) : op(mdl::compile_hamiltonian(h).op), basis(mdl::physical_basis(h)) {}
};

cd::Options options_for(std::size_t states, Sector sector) {
    auto options = exact_options(states, std::move(sector));
    return options;
}

cd::Evolution real_time(double dt, std::size_t steps, std::size_t states) {
    cd::Evolution evolution;
    evolution.time_step = dt;
    evolution.steps = steps;
    evolution.states = states;
    return evolution;
}

// The exact Trotter product of the engine's schedule, one step per call, a
// drive's coefficient being that of the step.
struct Trotter {
    Reference const& reference;
    Dense half, full;
    bool imaginary;
    double dt;

    Trotter(Reference const& r, mdl::CompiledOperator const& h, double time_step, bool imag = false)
        : reference(r), imaginary(imag), dt(time_step) {
        set(h);
    }

    void set(mdl::CompiledOperator const& h) {
        auto const [even, odd] = even_odd(h);
        Complex const unit = imaginary ? Complex(-1.) : Complex(0., -1.);
        half = exponential(reference.matrix(even), unit * (0.5 * dt));
        full = exponential(reference.matrix(odd), unit * dt);
    }

    Vector step(Vector const& psi) const {
        Vector out = half * (full * (half * psi));
        if (imaginary) out /= out.norm();
        return out;
    }
};

// min over phases of |a - e^{i phi} b|, without the cancellation of the
// closed form at small distances.
double distance_up_to_phase(Vector const& a, Vector const& b) {
    Complex const overlap = b.dot(a); // conj(b) . a
    Complex const phase = std::abs(overlap) > 0. ? overlap / std::abs(overlap) : Complex(1.);
    return (a - phase * b).norm();
}

// Evolve with every step measured and the states kept, and compare each step
// with the Trotter product (to `trotter_tolerance`) and the exact evolution
// (to `exact_tolerance`).
void expect_trotter_and_exact(Model const& initial, Model const& final_model, Sector const& sector,
                              cd::Measurements const& requests, cd::Evolution evolution, std::size_t prepare_states,
                              double trotter_tolerance, double exact_tolerance) {
    evolution.keep_mps = true;
    auto const result = cd::evolve(initial.op, final_model.op, initial.basis, evolution,
                                   options_for(prepare_states, sector), requests);
    Reference const reference(initial.op, initial.basis, sector);
    Trotter const trotter(reference, final_model.op, evolution.time_step);
    Dense const h = reference.matrix(final_model.op);
    ASSERT_EQ(result.steps.size(), evolution.steps + 1);
    ASSERT_EQ(result.mps.size(), evolution.steps + 1);
    ASSERT_EQ(result.truncation_errors.size(), evolution.steps + 1);
    EXPECT_NEAR(result.initial.energies.at(0), reference.states.eigenvalues.at(0), 1e-9);
    Vector psi = reference.initial;
    Vector const start = reference.vector(result.mps[0]);
    for (std::size_t k = 0; k <= evolution.steps; ++k) {
        SCOPED_TRACE("step " + std::to_string(k));
        if (k > 0) psi = trotter.step(psi);
        Vector const exact = exponential(h, Complex(0., -evolution.time_step * static_cast<double>(k))) * reference.initial;
        EXPECT_EQ(result.steps[k], k);
        EXPECT_DOUBLE_EQ(result.times[k], evolution.time_step * static_cast<double>(k));
        EXPECT_LT(result.truncation_errors[k], 1e-12);
        Vector const state = reference.vector(result.mps[k]);
        EXPECT_NEAR(state.norm(), 1., 1e-10);
        EXPECT_LT(distance_up_to_phase(psi, state), trotter_tolerance);
        EXPECT_LT(distance_up_to_phase(exact, state), exact_tolerance);
        // The Loschmidt amplitude <psi(0)|psi(t)> is phase-free.
        Complex const amplitude = start.dot(state), expected = reference.initial.dot(psi);
        EXPECT_LT(std::abs(amplitude - expected), trotter_tolerance);
        EXPECT_NEAR(result.energies[k], reference.value(psi, final_model.op), trotter_tolerance);
        for (auto const& [name, op] : requests.terms) {
            auto const values = reference.values(psi, op);
            auto const& rows = result.measurements.at(name).values;
            ASSERT_EQ(rows.size(), evolution.steps + 1);
            ASSERT_EQ(rows[k].size(), values.size());
            for (std::size_t c = 0; c < values.size(); ++c)
                EXPECT_NEAR(rows[k][c], values[c], trotter_tolerance) << name << " column " << c;
        }
    }
}

mdl::Hamiltonian chain_xxz(int length, double jxy, double jz, double h = 0.) {
    return mdl::xxz(lat::chain(length, 1.0, lat::Boundary::open), jxy, jz, h);
}

// Spinless fermions hopping with the phase exp(i phi) and a weak interaction.
mdl::Hamiltonian complex_hopping(int length, double phi) {
    auto const site = mdl::spinless_fermion();
    mdl::Hamiltonian model(lat::chain(length, 1.0, lat::Boundary::open), site, mdl::FermionConvention::strict, {"N"});
    model.add_bond(-std::polar(1.0, phi), [site](int i, int j) {
        return mdl::TermSum::located(site, "c_dag", i) * mdl::TermSum::located(site, "c", j);
    }, 0, true);
    model.add_bond(0.7, [site](int i, int j) {
        return mdl::TermSum::located(site, "n", i) * mdl::TermSum::located(site, "n", j);
    }, 0, false);
    return model;
}

} // namespace

TEST(ClassicDmrgEvolution, VanishingFinalHamiltonianKeepsTheState) {
    // Identity gates: the wave-function transformations alone carry the state
    // through several sweeps, measurements and the exported states unchanged.
    Model const initial(chain_xxz(8, 1.0, 0.5));
    mdl::CompiledOperator const nothing{initial.op.site_dims, {}};
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(chain_xxz(8, 1.0, 0.5), "Sz");
    requests.terms["SzSz"] = cd::pair_terms(chain_xxz(8, 1.0, 0.5), "Sz", "Sz");
    auto evolution = real_time(0.1, 4, 64);
    evolution.keep_mps = true;
    auto const result = cd::evolve(initial.op, nothing, initial.basis, evolution, options_for(64, {{"Sz", 0.}}), requests);
    Reference const reference(initial.op, initial.basis, {{"Sz", 0.}});
    ASSERT_EQ(result.mps.size(), 5u);
    Vector const start = reference.vector(result.mps[0]);
    EXPECT_NEAR(std::abs(reference.initial.dot(start)), 1., 1e-10);
    for (std::size_t k = 0; k < result.mps.size(); ++k) {
        SCOPED_TRACE(k);
        Vector const state = reference.vector(result.mps[k]);
        EXPECT_LT(std::abs(start.dot(state) - 1.), 1e-10);
        EXPECT_NEAR(result.energies[k], 0., 1e-12);
        for (auto const& [name, rows] : result.measurements)
            for (std::size_t c = 0; c < rows.values[k].size(); ++c)
                EXPECT_NEAR(rows.values[k][c], rows.values[0][c], 1e-10) << name << " " << c;
    }
    auto const expected = reference.values(reference.initial, requests.terms.at("SzSz"));
    for (std::size_t c = 0; c < expected.size(); ++c)
        EXPECT_NEAR(result.measurements.at("SzSz").values[0][c], expected[c], 1e-10) << c;
}

TEST(ClassicDmrgEvolution, EigenstateOnlyGainsAPhase) {
    Model const model(chain_xxz(8, 1.0, 0.7));
    auto evolution = real_time(0.1, 5, 64);
    evolution.keep_mps = true;
    auto const result = cd::evolve(model.op, model.op, model.basis, evolution, options_for(64, {{"Sz", 0.}}));
    Reference const reference(model.op, model.basis, {{"Sz", 0.}});
    double const e0 = reference.states.eigenvalues[0];
    Vector const start = reference.vector(result.mps[0]);
    for (std::size_t k = 0; k < result.mps.size(); ++k) {
        SCOPED_TRACE(k);
        // The Trotter error of an eigenstate is the only change: second order in dt.
        Vector const state = reference.vector(result.mps[k]);
        EXPECT_GT(std::abs(start.dot(state)), 1. - 1e-4);
        EXPECT_NEAR(result.energies[k], e0, 1e-3);
    }
}

TEST(ClassicDmrgEvolution, XxzQuenchFollowsTheTrotterProduct) {
    // XX ground state evolved under the XXZ chain at Jz = 2, in a field.
    Model const initial(chain_xxz(8, 1.0, 0.0));
    Model const final_model(chain_xxz(8, 1.0, 2.0, 0.3));
    auto const lattice_model = chain_xxz(8, 1.0, 0.0);
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(lattice_model, "Sz");
    requests.terms["SzSz"] = cd::pair_terms(lattice_model, "Sz", "Sz");
    requests.terms["SpSm"] = cd::pair_terms(lattice_model, "Sp", "Sm");
    expect_trotter_and_exact(initial, final_model, {{"Sz", 0.}}, requests, real_time(0.05, 20, 64), 64, 1e-9, 2e-2);
}

TEST(ClassicDmrgEvolution, OddChainAndFourSites) {
    // An odd length (the half sweeps end on different bonds) and the shortest
    // chain, whose one configuration holds every bond.
    for (int length : {7, 4}) {
        SCOPED_TRACE(length);
        Sector const sector{{"Sz", length % 2 ? 0.5 : 0.}};
        auto const lattice_model = chain_xxz(length, 1.0, 0.0);
        cd::Measurements requests;
        requests.terms["Sz"] = cd::site_terms(lattice_model, "Sz");
        expect_trotter_and_exact(Model(chain_xxz(length, 1.0, 0.0)), Model(chain_xxz(length, 0.8, 1.5, -0.2)), sector,
                                 requests, real_time(0.05, 8, 64), 64, 1e-9, 1e-2);
    }
}

TEST(ClassicDmrgEvolution, XxDomainWallMeltsAsTheExactEvolution) {
    // |up..up down..down>, the ground state of opposite fields on the two halves, under the XX chain.
    int const length = 10;
    auto const xx = chain_xxz(length, 1.0, 0.0);
    Model const final_model(xx);
    auto field = cd::site_terms(xx, "Sz");
    for (auto& term : field.terms) term.coefficient *= term.factors[0].site < length / 2 ? -1. : 1.;
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(xx, "Sz");
    Model initial(xx);
    initial.op = field;
    expect_trotter_and_exact(initial, final_model, {{"Sz", 0.}}, requests, real_time(0.05, 30, 32), 32, 1e-9, 2e-2);
}

TEST(ClassicDmrgEvolution, HubbardQuenchKeepsTheFermionSigns) {
    auto const lattice = lat::chain(6, 1.0, lat::Boundary::open);
    auto const free = mdl::hubbard(lattice, 1.0, 0.0);
    Model const initial(free);
    Model const final_model(mdl::hubbard(lattice, 1.0, 4.0, 0.2));
    cd::Measurements requests;
    requests.terms["n"] = cd::site_terms(free, "n");
    requests.terms["double"] = cd::site_terms(free, "n_updn");
    requests.terms["hop_up"] = cd::pair_terms(free, "c_dag_up", "c_up");
    expect_trotter_and_exact(initial, final_model, {{"N", 6.}, {"Sz", 0.}}, requests, real_time(0.05, 12, 64),
                             64, 1e-9, 2e-2);
}

TEST(ClassicDmrgEvolution, SpinlessFermionsWithInteraction) {
    auto const lattice = lat::chain(8, 1.0, lat::Boundary::open);
    auto const free = mdl::spinless_fermion_model(lattice, 1.0, 0.0, 0.0);
    // A potential step prepares an inhomogeneous state; V = 1.5 evolves it.
    auto potential = cd::site_terms(free, "n");
    for (auto& term : potential.terms) term.coefficient *= term.factors[0].site < 4 ? -1. : 1.;
    Model initial(free);
    initial.op = sum(initial.op, potential);
    Model const final_model(mdl::spinless_fermion_model(lattice, 1.0, 0.0, 1.5));
    cd::Measurements requests;
    requests.terms["n"] = cd::site_terms(free, "n");
    requests.terms["hop"] = cd::pair_terms(free, "c_dag", "c");
    expect_trotter_and_exact(initial, final_model, {{"N", 4.}}, requests, real_time(0.05, 12, 16), 16, 1e-9, 2e-2);
}

TEST(ClassicDmrgEvolution, ComplexHoppingPhases) {
    // A current-carrying quench: the gates are complex beyond a global phase.
    auto const lattice = lat::chain(8, 1.0, lat::Boundary::open);
    auto const free = mdl::spinless_fermion_model(lattice, 1.0, 0.0, 0.0);
    auto potential = cd::site_terms(free, "n");
    for (auto& term : potential.terms) term.coefficient *= 0.3 * term.factors[0].site;
    Model initial(free);
    initial.op = sum(initial.op, potential);
    cd::Measurements requests;
    requests.terms["n"] = cd::site_terms(free, "n");
    requests.terms["hop"] = cd::pair_terms(free, "c_dag", "c");
    expect_trotter_and_exact(initial, Model(complex_hopping(8, 0.4)), {{"N", 4.}}, requests, real_time(0.05, 12, 16),
                             16, 1e-9, 2e-2);
}

TEST(ClassicDmrgEvolution, TrotterErrorShrinksAsTheSquareOfTheStep) {
    Model const initial(chain_xxz(8, 1.0, 0.0));
    Model const final_model(chain_xxz(8, 1.0, 1.0, 0.2));
    Reference const reference(initial.op, initial.basis, {{"Sz", 0.}});
    Dense const h = reference.matrix(final_model.op);
    Vector const exact = exponential(h, Complex(0., -1.)) * reference.initial; // t = 1
    std::vector<double> errors;
    for (double dt : {0.2, 0.1, 0.05}) {
        auto evolution = real_time(dt, static_cast<std::size_t>(std::lround(1. / dt)), 64);
        evolution.keep_mps = true;
        evolution.measure_steps = {evolution.steps};
        auto const result = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options_for(64, {{"Sz", 0.}}));
        ASSERT_EQ(result.mps.size(), 1u);
        errors.push_back(distance_up_to_phase(exact, reference.vector(result.mps[0])));
    }
    for (std::size_t i = 1; i < errors.size(); ++i) {
        double const ratio = errors[i - 1] / errors[i];
        EXPECT_GT(ratio, 3.5) << i;
        EXPECT_LT(ratio, 4.5) << i;
    }
}

TEST(ClassicDmrgEvolution, RealTimeConservesTheEnergy) {
    Model const initial(chain_xxz(10, 1.0, 0.0));
    Model const final_model(chain_xxz(10, 1.0, 1.0));
    auto const result = cd::evolve(initial.op, final_model.op, initial.basis, real_time(0.02, 50, 32),
                                   options_for(32, {{"Sz", 0.}}));
    ASSERT_EQ(result.energies.size(), 51u);
    for (double e : result.energies) EXPECT_NEAR(e, result.energies[0], 2e-3);
    // The initial XX ground state is not an eigenstate of the final chain.
    EXPECT_GT(std::abs(result.energies[0] - result.initial.energies[0]), 0.1);
}

TEST(ClassicDmrgEvolution, ImaginaryTimeReachesTheGroundState) {
    Model const initial(chain_xxz(8, 1.0, 0.0));
    Model const final_model(chain_xxz(8, 1.0, 1.0, 0.1));
    double const e0 = ed_energies(chain_xxz(8, 1.0, 1.0, 0.1), {{"Sz", 0.}})[0];
    for (auto arithmetic : {cd::Arithmetic::real, cd::Arithmetic::complex}) {
        SCOPED_TRACE(static_cast<int>(arithmetic));
        auto evolution = real_time(0.05, 200, 64);
        evolution.imaginary = true;
        evolution.measure_steps = {0, 50, 200};
        evolution.keep_mps = true;
        auto options = options_for(64, {{"Sz", 0.}});
        options.arithmetic = arithmetic;
        auto const result = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options);
        EXPECT_EQ(result.initial.options.arithmetic, arithmetic);
        ASSERT_EQ(result.energies.size(), 3u);
        EXPECT_GT(result.energies[0], result.energies[1]);
        EXPECT_GT(result.energies[1], result.energies[2]);
        // tau = 10: exp(-tau gap) is negligible; the Trotter shift is O(dt^2).
        EXPECT_NEAR(result.energies[2], e0, 2e-3);
        EXPECT_EQ(result.mps.back().sites.front().complex.empty(), arithmetic == cd::Arithmetic::real);
        EXPECT_DOUBLE_EQ(result.times.back(), 10.);
    }
}

TEST(ClassicDmrgEvolution, ImaginaryTimeFollowsTheTrotterProduct) {
    Model const initial(chain_xxz(8, 1.0, 0.0));
    Model const final_model(chain_xxz(8, 1.0, 1.5, 0.2));
    auto evolution = real_time(0.1, 10, 64);
    evolution.imaginary = true;
    evolution.keep_mps = true;
    auto options = options_for(64, {{"Sz", 0.}});
    options.arithmetic = cd::Arithmetic::real;
    auto const result = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options);
    Reference const reference(initial.op, initial.basis, {{"Sz", 0.}});
    Trotter const trotter(reference, final_model.op, evolution.time_step, true);
    Vector psi = reference.initial;
    for (std::size_t k = 0; k <= evolution.steps; ++k) {
        SCOPED_TRACE(k);
        if (k > 0) psi = trotter.step(psi);
        EXPECT_LT(distance_up_to_phase(psi, reference.vector(result.mps[k])), 1e-9);
        EXPECT_NEAR(result.energies[k], reference.value(psi, final_model.op), 1e-9);
    }
}

TEST(ClassicDmrgEvolution, DrivesTakeTheirCoefficientPerStep) {
    // A staggered field switched on as sin(omega t), sampled at the step midpoints.
    int const length = 8;
    auto const xx = chain_xxz(length, 1.0, 0.0);
    Model const model(xx);
    auto staggered = cd::site_terms(xx, "Sz");
    for (auto& term : staggered.terms) term.coefficient *= term.factors[0].site % 2 ? -1. : 1.;
    auto evolution = real_time(0.05, 16, 64);
    evolution.keep_mps = true;
    std::vector<double> coefficients;
    for (std::size_t n = 0; n < evolution.steps; ++n)
        coefficients.push_back(0.8 * std::sin(3. * (static_cast<double>(n) + 0.5) * evolution.time_step));
    evolution.drives.push_back(cd::Drive{staggered, coefficients});
    Model initial(xx);
    initial.op = sum(initial.op, scaled(staggered, 0.3));
    auto const result = cd::evolve(initial.op, model.op, model.basis, evolution, options_for(64, {{"Sz", 0.}}));
    Reference const reference(initial.op, model.basis, {{"Sz", 0.}});
    Trotter trotter(reference, model.op, evolution.time_step);
    Vector psi = reference.initial;
    for (std::size_t k = 0; k <= evolution.steps; ++k) {
        SCOPED_TRACE(k);
        if (k > 0) {
            trotter.set(sum(model.op, scaled(staggered, coefficients[k - 1])));
            psi = trotter.step(psi);
        }
        EXPECT_LT(distance_up_to_phase(psi, reference.vector(result.mps[k])), 1e-9);
        // The energy is the final Hamiltonian's, without the drive.
        EXPECT_NEAR(result.energies[k], reference.value(psi, model.op), 1e-9);
    }
}

TEST(ClassicDmrgEvolution, MeasuredStepsSelectTheReportedStates) {
    Model const initial(chain_xxz(8, 1.0, 0.0));
    Model const final_model(chain_xxz(8, 1.0, 2.0));
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(chain_xxz(8, 1.0, 0.0), "Sz");
    auto evolution = real_time(0.05, 9, 64);
    auto const every = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options_for(64, {{"Sz", 0.}}),
                                  requests);
    evolution.measure_steps = {9, 3, 4, 3};
    auto const some = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options_for(64, {{"Sz", 0.}}),
                                 requests);
    EXPECT_EQ(some.steps, (std::vector<std::size_t>{3, 4, 9}));
    ASSERT_EQ(some.times.size(), 3u);
    EXPECT_DOUBLE_EQ(some.times[2], 0.45);
    EXPECT_EQ(some.truncation_errors.size(), 10u);
    EXPECT_TRUE(some.mps.empty());
    auto const& rows = some.measurements.at("Sz");
    ASSERT_EQ(rows.values.size(), 3u);
    ASSERT_EQ(rows.sites.size(), 8u);
    EXPECT_EQ(rows.sites[3], std::vector<int>{3});
    for (std::size_t k = 0; k < 3; ++k) {
        std::size_t const step = some.steps[k];
        // Merged and unmerged half steps give the same state when nothing is truncated.
        EXPECT_NEAR(some.energies[k], every.energies[step], 1e-9);
        for (std::size_t c = 0; c < 8; ++c)
            EXPECT_NEAR(rows.values[k][c], every.measurements.at("Sz").values[step][c], 1e-9) << step << " " << c;
    }
}

TEST(ClassicDmrgEvolution, TruncationIsReported) {
    // Few kept states on a quench that builds entanglement.
    Model const initial(chain_xxz(12, 1.0, 0.0));
    Model const final_model(chain_xxz(12, 1.0, 3.0));
    auto evolution = real_time(0.1, 10, 6);
    std::vector<cd::Iteration> iterations;
    auto const result = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options_for(32, {{"Sz", 0.}}),
                                   {}, [&](cd::Iteration const& row) { iterations.push_back(row); });
    ASSERT_EQ(result.truncation_errors.size(), 11u);
    double total = 0.;
    for (double w : result.truncation_errors) {
        EXPECT_GE(w, 0.);
        total += w;
    }
    EXPECT_GT(total, 1e-6);
    // The callback sees the preparation's iterations, then the evolution's without an energy.
    ASSERT_GT(iterations.size(), result.initial.iterations.size());
    for (std::size_t i = 0; i < result.initial.iterations.size(); ++i) EXPECT_TRUE(std::isfinite(iterations[i].energy));
    for (std::size_t i = result.initial.iterations.size(); i < iterations.size(); ++i)
        EXPECT_TRUE(std::isnan(iterations[i].energy));
    // A truncation error bound keeps fewer states than the cap where it can.
    evolution.states = 64;
    evolution.truncation_error = 1e-4;
    auto const bounded = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options_for(32, {{"Sz", 0.}}));
    for (double w : bounded.truncation_errors) EXPECT_LT(w, 1e-4 * 4 * 12);
}

TEST(ClassicDmrgEvolution, HamiltonianOverloadChecksTheBasis) {
    auto const a = chain_xxz(6, 1.0, 0.0);
    auto const b = chain_xxz(6, 1.0, 1.0);
    auto const result = cd::evolve(a, b, real_time(0.1, 2, 32), options_for(32, {{"Sz", 0.}}));
    EXPECT_EQ(result.energies.size(), 3u);
    auto const spin_one = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0, 0.0, 1.0);
    EXPECT_THROW(cd::evolve(a, spin_one, real_time(0.1, 2, 32), options_for(32, {{"Sz", 0.}})), std::invalid_argument);
}

TEST(ClassicDmrgEvolution, RejectsWhatTheGatesCannotApply) {
    Model const initial(chain_xxz(8, 1.0, 0.0));
    auto const options = options_for(32, {{"Sz", 0.}});
    auto expect_message = [](std::function<void()> const& call, std::string const& fragment) {
        try {
            call();
            ADD_FAILURE() << "no exception; expected " << fragment;
        } catch (std::invalid_argument const& e) {
            EXPECT_NE(std::string(e.what()).find(fragment), std::string::npos) << e.what();
        }
    };
    // Next-nearest neighbours and a periodic chain's closing bond.
    Model const frustrated(mdl::heisenberg(lat::nnn_chain(8, 1.0, lat::Boundary::open), 1.0));
    expect_message([&] { cd::evolve(initial.op, frustrated.op, initial.basis, real_time(0.1, 1, 32), options); },
                   "the final Hamiltonian couples sites 0 and 2");
    Model const ring(mdl::xxz(lat::chain(8, 1.0, lat::Boundary::periodic), 1.0, 1.0));
    expect_message([&] { cd::evolve(initial.op, ring.op, initial.basis, real_time(0.1, 1, 32), options); },
                   "nearest-neighbour bond gates in chain order");
    // A drive that leaves the sector, or with the wrong number of coefficients.
    auto evolution = real_time(0.1, 2, 32);
    evolution.drives.push_back(cd::Drive{cd::site_terms(chain_xxz(8, 1.0, 0.0), "Sx"), {1., 1.}});
    expect_message([&] { cd::evolve(initial.op, initial.op, initial.basis, evolution, options); },
                   "drive 0 does not conserve");
    evolution.drives[0] = cd::Drive{cd::site_terms(chain_xxz(8, 1.0, 0.0), "Sz"), {1.}};
    expect_message([&] { cd::evolve(initial.op, initial.op, initial.basis, evolution, options); },
                   "one coefficient per step");
    // The controls of the evolution.
    auto bad = real_time(0., 2, 32);
    expect_message([&] { cd::evolve(initial.op, initial.op, initial.basis, bad, options); }, "time_step");
    bad = real_time(0.1, 2, 32);
    bad.measure_steps = {3};
    expect_message([&] { cd::evolve(initial.op, initial.op, initial.basis, bad, options); }, "after the last step");
    bad = real_time(0.1, 2, 32);
    bad.truncation_error = 1.5;
    expect_message([&] { cd::evolve(initial.op, initial.op, initial.basis, bad, options); }, "truncation_error");
    auto two = options;
    two.num_eigenvalues = 2;
    expect_message([&] { cd::evolve(initial.op, initial.op, initial.basis, real_time(0.1, 1, 32), two); },
                   "num_eigenvalues");
    auto kept = options;
    kept.keep_mps = true;
    expect_message([&] { cd::evolve(initial.op, initial.op, initial.basis, real_time(0.1, 1, 32), kept); },
                   "Evolution::keep_mps");
    // A non-Hermitian final Hamiltonian.
    auto raising = cd::site_terms(chain_xxz(8, 1.0, 0.0), "Sz");
    raising.terms[0].coefficient = Complex(0., 1.);
    expect_message([&] { cd::evolve(initial.op, raising, initial.basis, real_time(0.1, 1, 32), options); },
                   "final Hamiltonian is not Hermitian");
    // Real arithmetic needs a real Hamiltonian in imaginary time.
    auto imaginary = real_time(0.1, 1, 32);
    imaginary.imaginary = true;
    auto real = exact_options(32, {{"N", 4.}});
    real.arithmetic = cd::Arithmetic::real;
    Model const fermions(mdl::spinless_fermion_model(lat::chain(8, 1.0, lat::Boundary::open)));
    Model const twisted(complex_hopping(8, 0.4));
    expect_message([&] { cd::evolve(fermions.op, twisted.op, fermions.basis, imaginary, real); },
                   "complex arithmetic");
}

// ------------------------------------------------------- time-step targeting

namespace {

cd::Evolution targeting(double dt, std::size_t steps, std::size_t states, bool imaginary = false) {
    auto evolution = real_time(dt, steps, states);
    evolution.method = cd::EvolutionMethod::targeting;
    evolution.imaginary = imaginary;
    return evolution;
}

// The final Hamiltonian without its identity part, which neither method
// applies (a phase in real time).
mdl::CompiledOperator without_phase(mdl::CompiledOperator op) {
    op.terms.erase(std::remove_if(op.terms.begin(), op.terms.end(), [](auto const& t) { return t.factors.empty(); }),
                   op.terms.end());
    return op;
}

// Evolve by time-step targeting with every step measured and the states
// kept, and compare each step with exp(-i t H) psi(0) (or the normalized
// exp(-t H) psi(0)) on the ED basis. There is no time-step error, so
// `tolerance` bounds the Krylov and truncation errors. In real time every
// iteration after the first half sweep reports <H_final>, which is
// conserved; the first half sweep's report NaN.
void expect_exact_targeting(Model const& initial, Model const& final_model, Sector const& sector,
                            cd::Measurements const& requests, cd::Evolution evolution, cd::Options const& options,
                            double tolerance) {
    evolution.keep_mps = true;
    std::vector<double> energies;
    auto const result = cd::evolve(initial.op, final_model.op, initial.basis, evolution, options, requests,
                                   [&](cd::Iteration const& row) { energies.push_back(row.energy); });
    Reference const reference(initial.op, initial.basis, sector);
    Dense const h = reference.matrix(without_phase(final_model.op));
    ASSERT_EQ(result.steps.size(), evolution.steps + 1);
    ASSERT_EQ(result.mps.size(), evolution.steps + 1);
    EXPECT_NEAR(result.initial.energies.at(0), reference.states.eigenvalues.at(0), 1e-9);
    Vector const start = reference.vector(result.mps[0]);
    for (std::size_t k = 0; k <= evolution.steps; ++k) {
        SCOPED_TRACE("step " + std::to_string(k));
        double const t = evolution.time_step * static_cast<double>(k);
        Vector exact = exponential(h, evolution.imaginary ? Complex(-t) : Complex(0., -t)) * reference.initial;
        exact /= exact.norm();
        Vector const state = reference.vector(result.mps[k]);
        EXPECT_NEAR(state.norm(), 1., 1e-10);
        EXPECT_LT(distance_up_to_phase(exact, state), tolerance);
        if (!evolution.imaginary) {
            Complex const amplitude = start.dot(state), expected = reference.initial.dot(exact);
            EXPECT_LT(std::abs(amplitude - expected), tolerance);
        }
        EXPECT_NEAR(result.energies[k], reference.value(exact, final_model.op), tolerance);
        for (auto const& [name, op] : requests.terms) {
            auto const values = reference.values(exact, op);
            auto const& row = result.measurements.at(name).values.at(k);
            ASSERT_EQ(row.size(), values.size());
            for (std::size_t c = 0; c < values.size(); ++c)
                EXPECT_NEAR(row[c], values[c], tolerance) << name << " column " << c;
        }
    }
    // The preparation's iterations come first; the first half sweep of the
    // evolution (length - 3 iterations) builds the final Hamiltonian.
    std::size_t const half = static_cast<std::size_t>(reference.states.problem.site_dims().size()) - 3;
    std::size_t const first = result.initial.iterations.size();
    ASSERT_GT(energies.size(), first + half);
    double const conserved = reference.value(reference.initial, final_model.op);
    for (std::size_t i = first; i < energies.size(); ++i) {
        if (i < first + half) {
            EXPECT_TRUE(std::isnan(energies[i])) << "iteration " << i;
        } else {
            ASSERT_TRUE(std::isfinite(energies[i])) << "iteration " << i;
            if (!evolution.imaginary) EXPECT_NEAR(energies[i], conserved, 1e-9) << "iteration " << i;
        }
    }
}

} // namespace

TEST(ClassicDmrgTargeting, FieldsToTheXxChainAreTheExactEvolution) {
    // The initial Hamiltonian has no S+ or S-: every block operator of the
    // final one is built by the evolution's first two half sweeps.
    int const length = 8;
    auto const xx = chain_xxz(length, 1.0, 0.0);
    auto field = cd::site_terms(xx, "Sz");
    for (auto& term : field.terms) term.coefficient *= term.factors[0].site < length / 2 ? -1. : 1.;
    Model initial(xx);
    initial.op = field;
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(xx, "Sz");
    expect_exact_targeting(initial, Model(xx), {{"Sz", 0.}}, requests, targeting(0.1, 10, 64),
                           options_for(64, {{"Sz", 0.}}), 1e-9);
}

TEST(ClassicDmrgTargeting, LongRangeAndPeriodicFinalHamiltonians) {
    Model const initial(chain_xxz(8, 1.0, 0.5));
    cd::Measurements requests;
    requests.terms["SzSz"] = cd::pair_terms(chain_xxz(8, 1.0, 0.5), "Sz", "Sz");
    Model const frustrated(mdl::heisenberg(lat::nnn_chain(8, 1.0, lat::Boundary::open), 1.0));
    expect_exact_targeting(initial, frustrated, {{"Sz", 0.}}, requests, targeting(0.1, 8, 64),
                           options_for(64, {{"Sz", 0.}}), 1e-9);
    Model const ring(mdl::xxz(lat::chain(8, 1.0, lat::Boundary::periodic), 1.0, 2.0));
    expect_exact_targeting(initial, ring, {{"Sz", 0.}}, requests, targeting(0.1, 8, 64),
                           options_for(64, {{"Sz", 0.}}), 1e-9);
}

TEST(ClassicDmrgTargeting, HubbardQuenchKeepsTheFermionSigns) {
    auto const lattice = lat::chain(6, 1.0, lat::Boundary::open);
    auto const free = mdl::hubbard(lattice, 1.0, 0.0);
    cd::Measurements requests;
    requests.terms["n"] = cd::site_terms(free, "n");
    requests.terms["hop_up"] = cd::pair_terms(free, "c_dag_up", "c_up");
    Sector const sector{{"N", 6.}, {"Sz", 0.}};
    expect_exact_targeting(Model(free), Model(mdl::hubbard(lattice, 1.0, 4.0, 0.2)), sector, requests,
                           targeting(0.1, 6, 64), options_for(64, sector), 1e-9);
}

TEST(ClassicDmrgTargeting, ComplexHoppingAndFourSites) {
    auto const lattice = lat::chain(8, 1.0, lat::Boundary::open);
    auto const free = mdl::spinless_fermion_model(lattice, 1.0, 0.0, 0.0);
    auto potential = cd::site_terms(free, "n");
    for (auto& term : potential.terms) term.coefficient *= 0.3 * term.factors[0].site;
    Model initial(free);
    initial.op = sum(initial.op, potential);
    cd::Measurements requests;
    requests.terms["hop"] = cd::pair_terms(free, "c_dag", "c");
    expect_exact_targeting(initial, Model(complex_hopping(8, 0.4)), {{"N", 4.}}, requests, targeting(0.1, 6, 64),
                           options_for(64, {{"N", 4.}}), 1e-9);
    Model const four(chain_xxz(4, 1.0, 0.0));
    expect_exact_targeting(four, Model(chain_xxz(4, 1.0, 2.0, 0.3)), {{"Sz", 0.}}, {}, targeting(0.1, 6, 16),
                           options_for(16, {{"Sz", 0.}}), 1e-9);
}

TEST(ClassicDmrgTargeting, ImaginaryTimeIsTheNormalizedExactProjection) {
    Model const initial(chain_xxz(8, 1.0, 0.0));
    Model const final_model(chain_xxz(8, 1.0, 1.0));
    expect_exact_targeting(initial, final_model, {{"Sz", 0.}}, {}, targeting(0.2, 10, 64, true),
                           options_for(64, {{"Sz", 0.}}), 1e-9);
    // In real arithmetic as well.
    auto options = options_for(64, {{"Sz", 0.}});
    options.arithmetic = cd::Arithmetic::real;
    expect_exact_targeting(initial, final_model, {{"Sz", 0.}}, {}, targeting(0.2, 10, 64, true), options, 1e-9);
}

TEST(ClassicDmrgTargeting, LargeStepsCarryNoTrotterError) {
    // dt = 0.4: the Krylov step is exact where the Trotter product is not.
    Model const initial(chain_xxz(8, 1.0, 0.0));
    Model const final_model(chain_xxz(8, 1.0, 2.0, 0.3));
    expect_exact_targeting(initial, final_model, {{"Sz", 0.}}, {}, targeting(0.4, 5, 64),
                           options_for(64, {{"Sz", 0.}}), 1e-9);
    auto trotter = real_time(0.4, 5, 64);
    trotter.keep_mps = true;
    auto const result = cd::evolve(initial.op, final_model.op, initial.basis, trotter, options_for(64, {{"Sz", 0.}}));
    Reference const reference(initial.op, initial.basis, {{"Sz", 0.}});
    Vector const exact = exponential(reference.matrix(final_model.op), Complex(0., -2.0)) * reference.initial;
    EXPECT_GT(distance_up_to_phase(exact, reference.vector(result.mps.back())), 1e-2);
}

TEST(ClassicDmrgTargeting, TargetsCarryAProductStateIntoTheEvolvedBasis) {
    // A domain wall is a product state: its density matrices have rank one,
    // so a truncation that targeted the state alone would keep one state per
    // block, and the Krylov step could move spins only on the two explicit
    // sites. The intermediate states of the step give the blocks the
    // directions the wall melts into.
    int const length = 10;
    auto const xx = chain_xxz(length, 1.0, 0.0);
    auto field = cd::site_terms(xx, "Sz");
    for (auto& term : field.terms) term.coefficient *= term.factors[0].site < length / 2 ? -1. : 1.;
    Model initial(xx);
    initial.op = field;
    auto options = options_for(64, {{"Sz", 0.}});
    options.warmup_configuration.assign(static_cast<std::size_t>(length), 1);
    for (int i = 0; i < length / 2; ++i) options.warmup_configuration[static_cast<std::size_t>(i)] = 0;
    auto evolution = targeting(0.05, 20, 64);
    evolution.truncation_error = 1e-12;
    cd::Measurements requests;
    requests.terms["Sz"] = cd::site_terms(xx, "Sz");
    expect_exact_targeting(initial, Model(xx), {{"Sz", 0.}}, requests, evolution, options, 1e-4);
}

TEST(ClassicDmrgTargeting, RejectsDrivesAndEmptySchedules) {
    Model const initial(chain_xxz(8, 1.0, 0.0));
    auto const options = options_for(32, {{"Sz", 0.}});
    auto evolution = targeting(0.1, 2, 32);
    evolution.drives.push_back(cd::Drive{cd::site_terms(chain_xxz(8, 1.0, 0.0), "Sz"), {1., 1.}});
    EXPECT_THROW(cd::evolve(initial.op, initial.op, initial.basis, evolution, options), std::invalid_argument);
    evolution = targeting(0.1, 2, 32);
    evolution.half_sweeps = 0;
    EXPECT_THROW(cd::evolve(initial.op, initial.op, initial.basis, evolution, options), std::invalid_argument);
}
