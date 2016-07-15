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

    alps::params p;
    set_common_parameters(p);
    CT_=new DoubleCellClusterTransformer(p);
    srand48(0);
  }
  void set_common_parameters(alps::params &p){
    p["dca.L"] = 20;
    p["dca.LATTICE"] = "Betts2D-8A";
    p["dca.SITES"] = n_sites_;
    p["dca.LATTICE_LIBRARY"] = "../../libcluster/cluster.xml";
    p["dca.LATTICE_DUMP"] = "./cluster.h5";
    p["MU"] = mu_;
    p["t"] = 1.;
    p["tprime"] = 0.;
    p["t_flavor"] = 0.;
    p["BETA"] = beta_;
    p["FLAVORS"] = 2;
    p["NMATSUBARA"] = n_freq_;
    p["N"] = n_tau_;
    p["dca.AFM_SELFCONSISTENCY"]=true;
    p["dca.PHASE"]=false;
    p["dca.IMBALANCE"]=false;
    p["dca.G0_OMEGA_INPUT"]="/tmp/none";
    p["dca.G0_OMEGA_OUTPUT"]="/tmp/none";
    p["dca.INITIAL_G0"]="FREE";
    p["dca.H_INIT"]=h_init_;
    p["dca.H"]=h_;
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
    for(int z=0;z<momentum        .size();++z) momentum        [z]=CT_->momentum_mesh().points()[k  ][z];
    for(int z=0;z<momentum_plus_pi.size();++z) momentum_plus_pi[z]=CT_->momentum_mesh().points()[k+1][z];
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
TEST_F(DoubleUnitCellTransformTest, ComputeNonintClusterSolution){
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
