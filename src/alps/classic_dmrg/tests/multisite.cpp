// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Numerical lock of the private dmtk engine on terms the public solve()
// does not accept: long-range hopping and four-site fermion products,
// against an independent occupation-bit exact diagonalization. The
// molecular solver relies on these engine paths. Adapted from
// skilledwolf/alpsx (private repository), commit 95cf817d6,
// tests/dmrg/classic_multisite_fermion_lock.cpp: a command-line program
// there; here GoogleTest cases for L = 4, 6, 8, run both with and without
// stored block products (the molecular solver stores them).
#include "invocation.hpp"
#include <gtest/gtest.h>

#include <Eigen/Eigenvalues>

#include <bitset>
#include <map>
#include <sstream>
#include <vector>

namespace {

// Set bits of an occupation word (C++17 has no std::popcount).
int popcount(unsigned bits) { return static_cast<int>(std::bitset<32>(bits).count()); }

struct Factor { int site; bool create; };
struct Term { double coefficient; std::vector<Factor> factors; };

void lock(int L, bool store_products) {
    int const N = L / 2;
    std::vector<Term> terms;
    for (int i = 0; i < L; ++i) {
        terms.push_back({.13 * (i + 1) * (i + 2), {{i, true}, {i, false}}});
        for (int j = i + 1; j < L; ++j) {
            terms.push_back({-.7 / (j - i), {{i, true}, {j, false}}});
            terms.push_back({-.7 / (j - i), {{j, true}, {i, false}}});
        }
    }
    terms.push_back({.37, {{0, true}, {1, true}, {L - 2, false}, {L - 1, false}}});
    terms.push_back({.37, {{L - 1, true}, {L - 2, true}, {1, false}, {0, false}}});
    std::vector<unsigned> states;
    std::map<unsigned, int> index;
    for (unsigned bits = 0; bits < (1u << L); ++bits)
        if (popcount(bits) == N) {
            index[bits] = static_cast<int>(states.size());
            states.push_back(bits);
        }
    auto const size = static_cast<Eigen::Index>(states.size());
    auto dense = [&](Term const& t) {
        Eigen::MatrixXd m = Eigen::MatrixXd::Zero(size, size);
        for (int col = 0; col < static_cast<int>(states.size()); ++col) {
            unsigned bits = states[static_cast<std::size_t>(col)];
            double coefficient = t.coefficient;
            for (auto p = t.factors.rbegin(); p != t.factors.rend(); ++p) {
                bool const occupied = (bits >> p->site) & 1u;
                if (occupied == p->create) { coefficient = 0.; break; }
                if (popcount(bits & ((1u << p->site) - 1u)) % 2) coefficient *= -1;
                bits ^= 1u << p->site;
            }
            if (coefficient) m(index.at(bits), col) += coefficient;
        }
        return m;
    };
    Eigen::MatrixXd exact = Eigen::MatrixXd::Zero(size, size);
    for (auto const& t : terms) exact += dense(t);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ed(exact);

    std::ostringstream log;
    dmtk::InvocationGuard guard({}, &log);
    dmtk::QN::add_qn_index("N", true);
    dmtk::QN::set_qn_mask(1);
    dmtk::Basis basis(2);
    for (int i = 0; i < 2; ++i) {
        dmtk::QN q;
        q["N"] = i;
        basis[static_cast<std::size_t>(i)] = dmtk::State(static_cast<std::size_t>(i), q);
    }
    basis.reorder();
    dmtk::Block<double> block;
    block.resize(basis);
    block.set_lattice(dmtk::Lattice(1, dmtk::OBC));
    auto op = [&](bool create) {
        dmtk::BasicOp<double> result(create ? "cd" : "c", 0);
        result.dqn["N"] = create ? 1 : -1;
        result.set_fermion(true).set_diagonal(false);
        result.resize(basis);
        dmtk::Matrix<double> matrix(2, 2);
        matrix = 0.;
        matrix[create ? 1 : 0][create ? 0 : 1] = 1.;
        result = matrix;
        block.push_back(result);
        return result;
    };
    auto cd = op(true), c = op(false);
    dmtk::BasicOp<double> n("n", 0);
    n.resize(basis);
    dmtk::Matrix<double> nm(2, 2);
    nm = 0.;
    nm[1][1] = 1.;
    n = nm;
    block.push_back(n);
    dmtk::Lattice lattice(L, dmtk::OBC);
    dmtk::Hami<double> h(lattice);
    for (int i = 0; i < L; ++i) h.sites(static_cast<std::size_t>(i)) = &block;
    auto native = [&](Term const& term) {
        dmtk::Term<double> t;
        if (term.factors.size() == 2 && term.factors[0].site == term.factors[1].site) {
            auto located = n;
            located.set_site(term.factors[0].site);
            t *= located;
        } else {
            // dmtk stores the rightmost factor first; operator*= appends.
            for (auto f = term.factors.rbegin(); f != term.factors.rend(); ++f) {
                auto located = f->create ? cd : c;
                located.set_site(f->site);
                t *= located;
            }
        }
        t.coef() = term.coefficient;
        return t.reorder();
    };
    for (auto const& term : terms) h += native(term);
    dmtk::System<double> s(h, lattice, "multisite_lock");
    s.set_store_products(store_products).set_lanczos_tolerance(1.e-12);
    s.qnt["N"] = N;
    s.set_qn_mask(1).set_use_hc(false).set_grow_symmetric(false).set_verbose(0);
    dmtk::Matrix<size_t> schedule(2, 5);
    schedule = 64;
    s.start(5, schedule);
    s.run(64);
    for (auto const& t : terms) s.corr += dmtk::BasicOp<double>(native(t));
    if (store_products) s.set_full_sweep(true);
    s.final_sweep(64, dmtk::RIGHT2LEFT, 1, true);
    s.measure();
    EXPECT_NEAR(s.energy[0], ed.eigenvalues()[0], 1e-9) << "L=" << L;
    for (std::size_t i = 0; i < terms.size(); ++i) {
        double const expected = ed.eigenvectors().col(0).dot(dense(terms[i]) * ed.eigenvectors().col(0));
        EXPECT_NEAR(s.corr[i].value(), expected, 1e-8) << "L=" << L << " term " << i;
    }
}

} // namespace

TEST(ClassicDmrgEngine, FourSiteFermionProducts) {
    for (int L : {4, 6, 8}) lock(L, false);
}

TEST(ClassicDmrgEngine, FourSiteFermionProductsWithStoredProducts) {
    for (int L : {4, 6, 8}) lock(L, true);
}
