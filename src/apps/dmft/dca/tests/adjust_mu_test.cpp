// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// sc.ADJUST_MU (matsubara_k_sum_adjust_mu): the density of the k-sum G is
// -G(beta), from its Fourier transform with a high-frequency tail. In the AFM
// and charge-order frameworks that tail had c2 = c3 = 0, so the MU found was
// wrong by an amount falling like 1/NMATSUBARA: 3e-2 at NMATSUBARA = 128 for
// Betts2D-4A at U = 0 (#449). The tail now follows from the patch moments of
// H and the tail of Sigma, which Dyson's equation gives in both frameworks.

#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace {
alps::params parameters(bool afm, int nmatsubara) {
  alps::params p;
  define_dca_parameters(p);
  p["dca.LATTICE"] = "Betts2D-4A";
  p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
  p["dca.L"] = 32;
  p["dca.SITES"] = 4;
  p["dca.AFM_SELFCONSISTENCY"] = afm;
  p["MU"] = 0.;
  p["t"] = 1.;
  p["tprime"] = 0.;
  p["BETA"] = 10.;
  p["FLAVORS"] = 2;
  p["NMATSUBARA"] = nmatsubara;
  p["N"] = 200;
  p["U"] = 0.;
  p["dca.G0_OMEGA_INPUT"] = "none";
  p["dca.G0_OMEGA_OUTPUT"] = "none";
  return p;
}

// Sigma(iw) = s0 + s1/(iw) for every momentum and spin, with its tail. In the
// double cell s0 and s1 are 2x2 matrices in (K, K+Q).
cluster_matsubara_kspace_gf normal_sigma(const NormalStateClusterTransformer& T, double s0, double s1) {
  cluster_matsubara_kspace_gf S = T.make_matsubara_gf();
  cluster_matsubara_kspace_gf::tail_type t0(S.mesh2(), S.mesh3()), t1(S.mesh2(), S.mesh3());
  for (alps::gf::matsubara_index w(0); w < S.mesh1().extent(); ++w)
    for (alps::gf::momentum_index k(0); k < S.mesh2().extent(); ++k)
      for (alps::gf::index f(0); f < S.mesh3().extent(); ++f)
        S(w, k, f) = s0 + s1 / std::complex<double>(0., S.mesh1().points()[w()]);
  for (alps::gf::momentum_index k(0); k < S.mesh2().extent(); ++k)
    for (alps::gf::index f(0); f < S.mesh3().extent(); ++f) {
      t0(k, f) = s0;
      t1(k, f) = s1;
    }
  S.set_tail(0, t0);
  S.set_tail(1, t1);
  return S;
}

double_cell_matsubara_kspace_gf double_cell_sigma(const DoubleCellClusterTransformer& T, const double s0[2][2], const double s1[2][2]) {
  double_cell_matsubara_kspace_gf S = T.make_matsubara_gf();
  double_cell_matsubara_kspace_gf::tail_type t0(S.mesh2(), S.mesh3(), S.mesh4()), t1(S.mesh2(), S.mesh3(), S.mesh4());
  for (alps::gf::momentum_index k(0); k < S.mesh2().extent(); ++k)
    for (alps::gf::index j(0); j < 2; ++j)
      for (alps::gf::index spin(0); spin < S.mesh4().extent(); ++spin) {
        const int i = k() % 2;
        t0(k, j, spin) = s0[i][j()];
        t1(k, j, spin) = s1[i][j()];
        for (alps::gf::matsubara_index w(0); w < S.mesh1().extent(); ++w)
          S(w, k, j, spin) = s0[i][j()] + s1[i][j()] / std::complex<double>(0., S.mesh1().points()[w()]);
      }
  S.set_tail(0, t0);
  S.set_tail(1, t1);
  return S;
}

double normal_mu(int nmatsubara, double target, double s0, double s1) {
  NormalStateClusterTransformer T(parameters(false, nmatsubara));
  cluster_matsubara_kspace_gf G = T.make_matsubara_gf(), S = normal_sigma(T, s0, s1);
  return T.matsubara_k_sum_adjust_mu(G, S, target);
}

double double_cell_mu(int nmatsubara, double target, const double s0[2][2], const double s1[2][2]) {
  DoubleCellClusterTransformer T(parameters(true, nmatsubara));
  double_cell_matsubara_kspace_gf G = T.make_matsubara_gf(), S = double_cell_sigma(T, s0, s1);
  return T.matsubara_k_sum_adjust_mu(G, S, target);
}
}

// At U = 0 the AFM framework without a field is the same lattice as the
// normal state, and so is a self-energy that is diagonal in (K, K+Q) and the
// same for K and K+Q: both frameworks have to find the same MU, at the same
// number of frequencies. On master the AFM MU was off by 3e-2 here.
TEST(AdjustMu, DoubleCellFindsTheMuOfTheNormalState) {
  const double zero[2][2] = {{0, 0}, {0, 0}};
  const double s0[2][2] = {{0.7, 0}, {0, 0.7}}, s1[2][2] = {{1.3, 0}, {0, 1.3}};
  for (double target : {0.6, 0.95}) {
    EXPECT_NEAR(double_cell_mu(128, target, zero, zero), normal_mu(128, target, 0., 0.), 1e-6) << "target " << target;
    EXPECT_NEAR(double_cell_mu(128, target, s0, s1), normal_mu(128, target, 0.7, 1.3), 1e-6) << "target " << target;
  }
}

// With a constant self-energy and a 1/iw term, also off-diagonal in the
// double cell, the MU depends on the number of frequencies only through the
// truncation of the Fourier transform, a few 1e-5 at NMATSUBARA = 128.
TEST(AdjustMu, NearlyIndependentOfTheNumberOfFrequencies) {
  const double target = 0.8;
  EXPECT_NEAR(normal_mu(128, target, 0.7, 1.3), normal_mu(1024, target, 0.7, 1.3), 1e-4);
  const double s0[2][2] = {{0.7, 0.3}, {0.3, 0.5}}, s1[2][2] = {{1.3, 0.2}, {0.2, 0.9}};
  EXPECT_NEAR(double_cell_mu(128, target, s0, s1), double_cell_mu(1024, target, s0, s1), 1e-4);
}

TEST(AdjustMu, NeedsTheTailOfSigma) {
  NormalStateClusterTransformer T(parameters(false, 64));
  cluster_matsubara_kspace_gf G = T.make_matsubara_gf(), S = T.make_matsubara_gf();
  S.initialize();
  try {
    T.matsubara_k_sum_adjust_mu(G, S, 0.8);
    FAIL() << "a self-energy without tail was accepted";
  } catch (std::runtime_error const& e) {
    EXPECT_NE(std::string(e.what()).find("the impurity solver returned a G without the tail orders 1 to 3"),
              std::string::npos) << e.what();
  }
}
