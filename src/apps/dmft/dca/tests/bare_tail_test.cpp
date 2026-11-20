// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// The bare tails attached to G0 must be the high-frequency expansion
// G0 = c1/iw + c2/(iw)^2 + c3/(iw)^3 + ... of the G0 they are attached to;
// the Fourier transform to G0(tau) relies on them. Each test computes G0 in
// one framework (normal state, double unit cell, Nambu) and one way (isolated
// cluster, free lattice, Dyson equation with a self-energy) and compares the
// tails with G0 at the highest Matsubara frequency.

#include <complex>
#include <memory>
#include <string>

#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"

namespace {

const int n_sites = 8;
const double beta = 5.;
const int n_freq = 4000;  // highest frequency ~5000

alps::params common_parameters() {
  alps::params p;
  p["dca.L"] = 8;
  p["dca.LATTICE"] = "Betts2D-8A";
  p["dca.SITES"] = n_sites;
  p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
  p["MU"] = 0.3;
  p["t"] = 1.;
  p["tprime"] = -0.25;
  p["U"] = 4.;
  p["V"] = 0.;
  p["BETA"] = beta;
  p["FLAVORS"] = 2;
  p["NMATSUBARA"] = n_freq;
  p["N"] = 100;
  p["dca.AFM_SELFCONSISTENCY"] = false;
  p["dca.CO_SELFCONSISTENCY"] = false;
  p["dca.G0_OMEGA_INPUT"] = "/tmp/none";
  p["dca.G0_OMEGA_OUTPUT"] = "/tmp/none";
  p["dca.INITIAL_G0"] = "FREE";
  p["dca.H_INIT"] = 0.;
  p["dca.H"] = 0.;
  p["dca.MU_STAG_INIT"] = 0.;
  p["dca.MU_STAG"] = 0.;
  p["dca.ETA1_INIT"] = 0.;
  p["dca.ETA2_INIT"] = 0.;
  p["dca.ETA1"] = 0.;
  p["dca.ETA2"] = 0.;
  return p;
}

// A self-energy Sigma0 + Sigma1/iw that differs between all components.
template<class GF> GF self_energy(const GF &like) {
  GF sigma(like);
  const auto &data = sigma.data();
  std::complex<double> *values = sigma.data().data();
  const std::size_t per_frequency = data.num_elements()/n_freq;
  for (int w=0; w<n_freq; ++w) {
    const std::complex<double> iw(0., (2*w+1)*M_PI/beta);
    for (std::size_t j=0; j<per_frequency; ++j) {
      const double s0 = 0.5 + 0.1*(j % 7), s1 = 1.5 + 0.2*(j % 5);
      values[w*per_frequency + j] = s0 + s1/iw;
    }
  }
  return sigma;
}

// G0 - c1/iw - c2/(iw)^2 - c3/(iw)^3 is of order c4/w^4, with c4 = <H^3>
// <~ 6^3 here, so at w ~ 5000 the c2 check is left with ~1e-5 and the c3
// check with ~0.04. A wrong sign of a field in c2 is off by >= 0.3, and a
// missing c3 by about the bandwidth squared.
const double tol_c2 = 1.e-3;
const double tol_c3 = 0.3;

void expect_tail(const cluster_matsubara_kspace_gf &G0) {
  const alps::gf::matsubara_index w(n_freq-1);
  const std::complex<double> iw(0., G0.mesh1().points()[w()]);
  for (alps::gf::momentum_index k(0); k<G0.mesh2().extent(); ++k) {
    for (alps::gf::index f(0); f<G0.mesh3().extent(); ++f) {
      const double c1 = G0.tail(1)(k,f), c2 = G0.tail(2)(k,f), c3 = G0.tail(3)(k,f);
      const std::complex<double> g = G0(w,k,f);
      EXPECT_NEAR(std::abs(iw*iw*(g - c1/iw) - c2 - c3/iw), 0., tol_c2) << "c2 at k=" << k() << " f=" << f();
      EXPECT_NEAR(std::abs(iw*iw*iw*(g - c1/iw - c2/(iw*iw)) - c3), 0., tol_c3) << "c3 at k=" << k() << " f=" << f();
    }
  }
}

void expect_tail(const cluster_nambu_matsubara_kspace_gf &G0) {
  const alps::gf::matsubara_index w(n_freq-1);
  const std::complex<double> iw(0., G0.mesh1().points()[w()]);
  for (alps::gf::momentum_index k(0); k<G0.mesh2().extent(); ++k) {
    for (alps::gf::index i(0); i<2; ++i) {
      for (alps::gf::index j(0); j<2; ++j) {
        const double c1 = G0.tail(1)(k,i,j), c2 = G0.tail(2)(k,i,j), c3 = G0.tail(3)(k,i,j);
        const std::complex<double> g = G0(w,k,i,j);
        EXPECT_NEAR(std::abs(iw*iw*(g - c1/iw) - c2 - c3/iw), 0., tol_c2)
            << "c2 at k=" << k() << " (" << i() << "," << j() << ")";
        EXPECT_NEAR(std::abs(iw*iw*iw*(g - c1/iw - c2/(iw*iw)) - c3), 0., tol_c3)
            << "c3 at k=" << k() << " (" << i() << "," << j() << ")";
      }
    }
  }
}

}  // namespace

class NormalStateBareTail : public ::testing::Test {
protected:
  void SetUp() override {
    CT_.reset(new NormalStateClusterTransformer(common_parameters()));
  }
  cluster_matsubara_kspace_gf make_gf() const {
    return cluster_matsubara_kspace_gf(alps::gf::omega_k_sigma_gf(
        alps::gf::matsubara_positive_mesh(beta, n_freq), CT_->momentum_mesh(), alps::gf::index_mesh(2)));
  }
  std::unique_ptr<NormalStateClusterTransformer> CT_;
};
TEST_F(NormalStateBareTail, IsolatedCluster) {
  cluster_matsubara_kspace_gf G0 = make_gf();
  CT_->compute_isolated_cluster_G0(G0);
  expect_tail(G0);
}
TEST_F(NormalStateBareTail, FreeLattice) {
  cluster_matsubara_kspace_gf G0 = make_gf();
  CT_->compute_noninteracting_G0(G0);
  expect_tail(G0);
}
TEST_F(NormalStateBareTail, Dyson) {
  cluster_matsubara_kspace_gf G = make_gf(), G0 = make_gf();
  const cluster_matsubara_kspace_gf sigma = self_energy(make_gf());
  CT_->matsubara_k_sum(G, sigma);
  CT_->dyson_get_green0(G0, sigma, G);
  expect_tail(G0);
}

// AFM: a staggered field h couples k and k+Q; CO: a staggered chemical potential.
class DoubleCellBareTail : public ::testing::TestWithParam<bool> {
protected:
  void SetUp() override {
    alps::params p = common_parameters();
    const bool afm = GetParam();
    p["dca.AFM_SELFCONSISTENCY"] = afm;
    p["dca.CO_SELFCONSISTENCY"] = !afm;
    if (afm) { p["dca.H_INIT"] = 0.4; p["dca.H"] = 0.25; }
    else     { p["dca.MU_STAG_INIT"] = 0.35; p["dca.MU_STAG"] = 0.15; }
    CT_.reset(new DoubleCellClusterTransformer(p));
  }
  double_cell_matsubara_kspace_gf make_gf() const {
    return double_cell_matsubara_kspace_gf(alps::gf::omega_k_sigma1_sigma2_gf(
        alps::gf::matsubara_positive_mesh(beta, n_freq), CT_->momentum_mesh(),
        alps::gf::index_mesh(2), alps::gf::index_mesh(2)));
  }
  std::unique_ptr<DoubleCellClusterTransformer> CT_;
};
TEST_P(DoubleCellBareTail, IsolatedCluster) {
  double_cell_matsubara_kspace_gf G0 = make_gf();
  CT_->compute_isolated_cluster_G0(G0);
  expect_tail(G0);
}
TEST_P(DoubleCellBareTail, FreeLattice) {
  double_cell_matsubara_kspace_gf G0 = make_gf();
  CT_->compute_noninteracting_G0(G0);
  expect_tail(G0);
}
TEST_P(DoubleCellBareTail, Dyson) {
  double_cell_matsubara_kspace_gf G = make_gf(), G0 = make_gf();
  const double_cell_matsubara_kspace_gf sigma = self_energy(make_gf());
  CT_->matsubara_k_sum(G, sigma);
  CT_->dyson_get_green0(G0, sigma, G);
  expect_tail(G0);
}
INSTANTIATE_TEST_SUITE_P(Order, DoubleCellBareTail, ::testing::Values(true, false),
                         [](const ::testing::TestParamInfo<bool> &info) { return info.param ? "AFM" : "CO"; });

// A Zeeman field h and a real d-wave pairing field eta1.
class NambuBareTail : public ::testing::Test {
protected:
  void SetUp() override {
    alps::params p = common_parameters();
    p["sc.NAMBU"] = true;
    p["dca.H_INIT"] = 0.3;
    p["dca.H"] = 0.2;
    p["dca.ETA1_INIT"] = 0.4;
    p["dca.ETA1"] = 0.15;
    CT_.reset(new NambuClusterTransformer(p));
  }
  cluster_nambu_matsubara_kspace_gf make_gf() const {
    return cluster_nambu_matsubara_kspace_gf(alps::gf::omega_k_sigma1_sigma2_gf(
        alps::gf::matsubara_positive_mesh(beta, n_freq), CT_->momentum_mesh(),
        alps::gf::index_mesh(2), alps::gf::index_mesh(2)));
  }
  std::unique_ptr<NambuClusterTransformer> CT_;
};
TEST_F(NambuBareTail, IsolatedCluster) {
  cluster_nambu_matsubara_kspace_gf G0 = make_gf();
  CT_->compute_isolated_cluster_G0(G0);
  expect_tail(G0);
}
TEST_F(NambuBareTail, FreeLattice) {
  cluster_nambu_matsubara_kspace_gf G0 = make_gf();
  CT_->compute_noninteracting_G0(G0);
  expect_tail(G0);
}
TEST_F(NambuBareTail, Dyson) {
  cluster_nambu_matsubara_kspace_gf G = make_gf(), G0 = make_gf();
  const cluster_nambu_matsubara_kspace_gf sigma = self_energy(make_gf());
  CT_->matsubara_k_sum(G, sigma);
  CT_->dyson_get_green0(G0, sigma, G);
  expect_tail(G0);
}
