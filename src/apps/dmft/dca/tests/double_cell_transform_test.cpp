// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"
#include<iostream>

class DoubleUnitCellTransformTest : public ::testing::Test{
protected:
  virtual void SetUp() {
    n_sites_=8;
    beta_=40;
    mu_=0.1;
    n_freq_=20;
    n_tau_=40;
    n_flavors_=2;
    h_init_=0.2;
    h_=0.1;
    mu_stag_init_=0.;
    mu_stag_=0.;

    alps::params p;
    set_common_parameters(p);
    CT_=new DoubleCellClusterTransformer(p);
    srand48(0);
  }
  void set_common_parameters(alps::params &p){
    p["dca.L"] = 20;
    p["dca.LATTICE"] = "Betts2D-8A";
    p["dca.SITES"] = n_sites_;
    p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
    p["MU"] = mu_;
    p["t"] = 1.;
    p["U"] = 3.1415;
    p["V"] = 1.2;
    p["tprime"] = 0.;
    p["BETA"] = beta_;
    p["FLAVORS"] = 2;
    p["NMATSUBARA"] = n_freq_;
    p["N"] = n_tau_;
    p["dca.AFM_SELFCONSISTENCY"]=true;
    p["dca.CO_SELFCONSISTENCY"]=false;
    p["dca.G0_OMEGA_INPUT"]="/tmp/none";
    p["dca.G0_OMEGA_OUTPUT"]="/tmp/none";
    p["dca.INITIAL_G0"]="FREE";
    p["dca.H_INIT"]=h_init_;
    p["dca.H"]=h_;
    p["dca.MU_STAG_INIT"]=mu_stag_init_;
    p["dca.MU_STAG"]=mu_stag_;
  }
  int L_;
  DoubleCellClusterTransformer *CT_;
  int n_freq_;
  int n_flavors_;
  int n_sites_;
  int n_tau_;
  double beta_;
  double mu_;
  double h_init_;
  double h_;
  double mu_stag_init_;
  double mu_stag_;

};
TEST_F(DoubleUnitCellTransformTest, MakeADoubleUnitCellTransformer){
;
}
TEST_F(DoubleUnitCellTransformTest, ReadNonexistentInputFile){
  double_cell_matsubara_kspace_gf G0_init(alps::gf::omega_k_sigma1_sigma2_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));
  EXPECT_EQ(CT_->read_initial_G0(G0_init), false);
}
TEST_F(DoubleUnitCellTransformTest, ComputeIsolatedClusterSolution){
  double_cell_matsubara_kspace_gf G0_init(alps::gf::omega_k_sigma1_sigma2_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));
  CT_->compute_isolated_cluster_G0(G0_init);
  for(int k=0;k<n_sites_;k+=2){
    std::vector<double> momentum(CT_->dimension());
    std::vector<double> momentum_plus_pi(CT_->dimension());
    for(int z=0;z<momentum        .size();++z) momentum        [z]=CT_->momentum_mesh().points()(k  ,z);
    for(int z=0;z<momentum_plus_pi.size();++z) momentum_plus_pi[z]=CT_->momentum_mesh().points()(k+1,z);
    for(int i=0;i<n_freq_;++i){
      for(int spin=0;spin<2;++spin){
        double spin_factor=spin==0?1.:-1.;
        //this implements Eq. 119
        two_by_two G(std::complex<double>(mu_-CT_->epsilon(momentum),(2*i+1)*M_PI/beta_),
            spin_factor*h_init_/2.,
            spin_factor*h_init_/2.,
            std::complex<double>(mu_-CT_->epsilon(momentum_plus_pi),(2*i+1)*M_PI/beta_));
        G.invert();
        EXPECT_NEAR(G(0,0).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(0), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G(0,1).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(1), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G(1,0).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(0), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G(1,1).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(1), alps::gf::index(spin)).real(), 1.e-10);
      }
    }
  }
}

TEST(TwoByTwo, AccessAndInversion){
  two_by_two G(1,2,3,4);
  EXPECT_NEAR(G(0,0).real(), 1, 1.e-10);
  EXPECT_NEAR(G(0,1).real(), 2, 1.e-10);
  EXPECT_NEAR(G(1,0).real(), 3, 1.e-10);
  EXPECT_NEAR(G(1,1).real(), 4, 1.e-10);
  G.invert();
  EXPECT_NEAR(G(0,0).real(), -2, 1.e-10);
  EXPECT_NEAR(G(0,1).real(), 1, 1.e-10);
  EXPECT_NEAR(G(1,0).real(), 1.5, 1.e-10);
  EXPECT_NEAR(G(1,1).real(), -0.5, 1.e-10);
}
TEST(TwoByTwo, Multiply){
  two_by_two G(1,2,3,4);
  two_by_two F(5,6,7,8);
  two_by_two H=G*F;
  EXPECT_NEAR(H(0,0).real(), 19, 1.e-10);
  EXPECT_NEAR(H(0,1).real(), 22, 1.e-10);
  EXPECT_NEAR(H(1,0).real(), 43, 1.e-10);
  EXPECT_NEAR(H(1,1).real(), 50, 1.e-10);
  two_by_two Z=G-H;
  EXPECT_NEAR(Z(0,0).real(), -18, 1.e-10);
  EXPECT_NEAR(Z(0,1).real(), -20, 1.e-10);
  EXPECT_NEAR(Z(1,0).real(), -40, 1.e-10);
  EXPECT_NEAR(Z(1,1).real(), -46, 1.e-10);
}
TEST(TwoByTwo, Scale){
  two_by_two G(1,2,3,4);
  two_by_two H=G*1.2;
  EXPECT_NEAR(H(0,0).real(), 1.2, 1.e-10);
  EXPECT_NEAR(H(0,1).real(), 2.4, 1.e-10);
  EXPECT_NEAR(H(1,0).real(), 3.6, 1.e-10);
  EXPECT_NEAR(H(1,1).real(), 4.8, 1.e-10);
}

TEST_F(DoubleUnitCellTransformTest, CompareIsolatedClusterSolutionToNormalStateForH0){
  alps::params p;
  set_common_parameters(p);
  p["dca.H_INIT"]=0;
  DoubleCellClusterTransformer CT1(p);
  NormalStateClusterTransformer CT2(p);

  double_cell_matsubara_kspace_gf G0_init(alps::gf::omega_k_sigma1_sigma2_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT1.momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));
  CT1.compute_isolated_cluster_G0(G0_init);

  cluster_matsubara_kspace_gf G0_init2(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT1.momentum_mesh(), alps::gf::index_mesh(2)));
  CT2.compute_isolated_cluster_G0(G0_init2);

  for(int k=0;k<n_sites_;k+=2){
    for(int i=0;i<n_freq_;++i){
      for(int spin=0;spin<2;++spin){
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(0), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).imag(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(0), alps::gf::index(spin)).imag(), 1.e-10);
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(spin)).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(1), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(spin)).imag(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(1), alps::gf::index(spin)).imag(), 1.e-10);
        EXPECT_NEAR(0., G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(1), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(0., G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(0), alps::gf::index(spin)).real(), 1.e-10);
      }
    }
  }
}
TEST_F(DoubleUnitCellTransformTest, CompareNonintClusterSolutionToNormalStateForH0){
  alps::params p;
  set_common_parameters(p);
  p["dca.H_INIT"]=0;
  DoubleCellClusterTransformer CT1(p);
  NormalStateClusterTransformer CT2(p);

  double_cell_matsubara_kspace_gf G0_init(alps::gf::omega_k_sigma1_sigma2_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT1.momentum_mesh(), alps::gf::index_mesh(2), alps::gf::index_mesh(2)));
  CT1.compute_noninteracting_G0(G0_init);

  cluster_matsubara_kspace_gf G0_init2(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT1.momentum_mesh(), alps::gf::index_mesh(2)));
  CT2.compute_noninteracting_G0(G0_init2);

  for(int k=0;k<n_sites_;k+=2){
    for(int i=0;i<n_freq_;++i){
      for(int spin=0;spin<2;++spin){
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(0), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).imag(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(0), alps::gf::index(spin)).imag(), 1.e-10);
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(spin)).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(1), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G0_init2(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(spin)).imag(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(1), alps::gf::index(spin)).imag(), 1.e-10);
        EXPECT_NEAR(0., G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(1), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(0., G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k+1), alps::gf::index(0), alps::gf::index(spin)).real(), 1.e-10);
      }
    }
  }
}

// The double unit cell needs a bipartite cluster: its blocks are (K, K+Q).
// Clusters that are not bipartite were accepted, and with an odd number of
// sites indexed past the momenta (Betts2D-1A crashed in the charge-order
// framework) (#451).
TEST(DoubleCellFramework, RefusesClustersThatAreNotBipartite){
  for (const char* framework : {"dca.AFM_SELFCONSISTENCY", "dca.CO_SELFCONSISTENCY"}) {
    struct cluster { const char* name; int sites; bool bipartite; };
    for (cluster c : {cluster{"Betts2D-4A", 4, true}, cluster{"Betts2D-8A", 8, true}, cluster{"Betts2D-16B", 16, true},
                      cluster{"Betts2D-1A", 1, false}, cluster{"Betts2D-9A", 9, false}, cluster{"Betts2D-10A", 10, false},
                      cluster{"Betts2D-13A", 13, false}}) {
      SCOPED_TRACE(std::string(framework) + " " + c.name);
      alps::params p;
      define_dca_parameters(p);
      p["dca.L"] = 4;
      p["dca.LATTICE"] = c.name;
      p["dca.SITES"] = c.sites;
      p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
      p["MU"] = 0.;
      p["t"] = 1.;
      p["U"] = 0.;
      p["tprime"] = 0.;
      p["BETA"] = 5.;
      p["FLAVORS"] = 2;
      p["NMATSUBARA"] = 8;
      p["N"] = 10;
      p[framework] = true;
      p["dca.G0_OMEGA_INPUT"] = "none";
      p["dca.G0_OMEGA_OUTPUT"] = "none";
      if (c.bipartite) {
        EXPECT_NO_THROW(DoubleCellClusterTransformer T(p));
      } else {
        try {
          DoubleCellClusterTransformer T(p);
          ADD_FAILURE() << "accepted";
        } catch (std::invalid_argument const& e) {
          const std::string what = e.what();
          EXPECT_NE(what.find(std::string(framework) + " needs a bipartite cluster"), std::string::npos) << what;
          EXPECT_NE(what.find("has " + std::to_string(c.sites) + " sites and is not bipartite"), std::string::npos) << what;
        }
      }
    }
  }
}

// symmetrize_k_space of the double unit cell was empty, so AFM and
// charge-order runs were not symmetrized over the cluster's point group
// (#453). It now averages over the classes of momentum pairs in a block.
namespace {
alps::params double_cell_parameters(const char* lattice, int sites, const char* framework, double h) {
  alps::params p;
  define_dca_parameters(p);
  p["dca.L"] = 8;
  p["dca.LATTICE"] = lattice;
  p["dca.SITES"] = sites;
  p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
  p["MU"] = 0.3;
  p["t"] = 1.;
  p["tprime"] = -0.2;
  p["U"] = 0.;
  p["BETA"] = 5.;
  p["FLAVORS"] = 2;
  p["NMATSUBARA"] = 6;
  p["N"] = 10;
  if (framework) p[framework] = true;
  if (framework && std::string(framework) == "dca.AFM_SELFCONSISTENCY") p["dca.H"] = h;
  if (framework && std::string(framework) == "dca.CO_SELFCONSISTENCY") p["dca.MU_STAG"] = h;
  p["dca.G0_OMEGA_INPUT"] = "none";
  p["dca.G0_OMEGA_OUTPUT"] = "none";
  return p;
}
std::complex<double> noise() { return std::complex<double>(drand48() - 0.5, drand48() - 0.5); }
}

// At zero field the AFM framework describes the paramagnet: with the same
// noise on the diagonal, its diagonal blocks are symmetrized as the normal
// state's G, and the off-diagonal blocks stay zero.
TEST(DoubleCellSymmetrization, ZeroFieldIsTheNormalState){
  for (auto c : {std::make_pair("Betts2D-8A", 8), std::make_pair("Betts2D-16B", 16)}) {
    SCOPED_TRACE(c.first);
    srand48(3);
    DoubleCellClusterTransformer D(double_cell_parameters(c.first, c.second, "dca.AFM_SELFCONSISTENCY", 0.));
    NormalStateClusterTransformer N(double_cell_parameters(c.first, c.second, nullptr, 0.));
    double_cell_matsubara_kspace_gf G = D.make_matsubara_gf();
    cluster_matsubara_kspace_gf Gn = N.make_matsubara_gf();
    G.initialize();
    for (alps::gf::matsubara_index w(0); w < G.mesh1().extent(); ++w)
      for (alps::gf::momentum_index k(0); k < G.mesh2().extent(); ++k)
        for (alps::gf::index s(0); s < 2; ++s) {
          const std::complex<double> x = noise();
          G(w, k, alps::gf::index(k() % 2), s) = x;
          Gn(w, k, s) = x;
        }
    D.symmetrize_k_space(G);
    N.symmetrize_k_space(Gn);
    for (alps::gf::matsubara_index w(0); w < G.mesh1().extent(); ++w)
      for (alps::gf::momentum_index k(0); k < G.mesh2().extent(); ++k)
        for (alps::gf::index s(0); s < 2; ++s) {
          EXPECT_NEAR(std::abs(G(w, k, alps::gf::index(k() % 2), s) - Gn(w, k, s)), 0., 1e-14);
          EXPECT_EQ(G(w, k, alps::gf::index(1 - k() % 2), s), std::complex<double>(0.));
        }
  }
}

// With a field and noise on every element, the symmetrized G is invariant
// under every symmetry of the cluster, G(gk1, gk2) = G(k1, k2), and a second
// symmetrization does not change it. In the charge-order framework also
// G(K, K+Q) = G(K+Q, K).
TEST(DoubleCellSymmetrization, InvariantAndIdempotent){
  for (const char* framework : {"dca.AFM_SELFCONSISTENCY", "dca.CO_SELFCONSISTENCY"}) {
    for (auto c : {std::make_pair("Betts2D-8A", 8), std::make_pair("Betts2D-16B", 16)}) {
      SCOPED_TRACE(std::string(framework) + " " + c.first);
      srand48(5);
      DoubleCellClusterTransformer D(double_cell_parameters(c.first, c.second, framework, 0.3));
      double_cell_matsubara_kspace_gf G = D.make_matsubara_gf();
      for (alps::gf::matsubara_index w(0); w < G.mesh1().extent(); ++w)
        for (alps::gf::momentum_index k(0); k < G.mesh2().extent(); ++k)
          for (alps::gf::index j(0); j < 2; ++j)
            for (alps::gf::index s(0); s < 2; ++s)
              G(w, k, j, s) = noise();
      D.symmetrize_k_space(G);
      double_cell_matsubara_kspace_gf G2 = G;
      D.symmetrize_k_space(G2);
      const auto& table = D.symmetry_table_k_space();
      EXPECT_GT(table.rows(), 1);
      auto element = [&](const double_cell_matsubara_kspace_gf& g, alps::gf::matsubara_index w, int k1, int k2, alps::gf::index s) {
        return g(w, alps::gf::momentum_index(k1), alps::gf::index(k2 % 2), s);
      };
      for (alps::gf::matsubara_index w(0); w < G.mesh1().extent(); ++w)
        for (int k1 = 0; k1 < c.second; ++k1)
          for (int k2 : {k1, k1 ^ 1})
            for (alps::gf::index s(0); s < 2; ++s) {
              EXPECT_NEAR(std::abs(element(G2, w, k1, k2, s) - element(G, w, k1, k2, s)), 0., 1e-14);
              for (int g = 0; g < table.rows(); ++g)
                EXPECT_NEAR(std::abs(element(G, w, table(g, k1), table(g, k2), s) - element(G, w, k1, k2, s)), 0., 1e-14)
                    << "symmetry " << g << ", pair " << k1 << " " << k2;
              if (std::string(framework) == "dca.CO_SELFCONSISTENCY")
                EXPECT_NEAR(std::abs(element(G, w, k2, k1, s) - element(G, w, k1, k2, s)), 0., 1e-14);
            }
    }
  }
}
