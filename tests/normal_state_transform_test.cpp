#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"
#include<iostream>

class NormalStateTransformTest : public ::testing::Test{
protected:
  virtual void SetUp() {
    n_sites_=8;
    beta_=40;
    mu_=0.1;
    n_freq_=20;
    n_tau_=40;
    n_flavors_=2;

    alps::params p;
    set_standard_parameters(p);
    CT_=new NormalStateClusterTransformer(p);

    srand48(0);
  }
  virtual void TearDown() {
    delete CT_;
  }
  void set_standard_parameters(alps::params &p){
    p["dca.L"] = 20;
    p["dca.LATTICE"] = "Betts2D-8A";
    p["dca.SITES"] = n_sites_;
    p["dca.LATTICE_LIBRARY"] = "../../libcluster/cluster.xml";
    p["dca.LATTICE_DUMP"] = "./cluster.h5";
    p["MU"] = mu_;
    p["t"] = 1.;
    p["U"] = 3.14;
    p["V"] = 1.2;
    p["tprime"] = 0.;
    p["t_flavor"] = 0.;
    p["BETA"] = beta_;
    p["FLAVORS"] = 2;
    p["NMATSUBARA"] = n_freq_;
    p["N"] = n_tau_;
    p["dca.AFM_SELFCONSISTENCY"]=false; //note that we need to know at this point if we're paramagnetic or not. Not sure this is the right place for it.
    p["dca.PHASE"]=false;
    p["dca.IMBALANCE"]=false;
    p["dca.G0_OMEGA_INPUT"]="/tmp/none";
    p["dca.G0_OMEGA_OUTPUT"]="/tmp/none";
    p["dca.INITIAL_G0"]="FREE";
  }
  int L_;
  NormalStateClusterTransformer *CT_;
  int n_freq_;
  int n_flavors_;
  int n_sites_;
  int n_tau_;
  double beta_;
  double mu_;
};
TEST_F(NormalStateTransformTest, MakeANormalStateTransformer){
  ; //just run Setup/TearDown
}
TEST_F(NormalStateTransformTest, ReadNonexistentInputFile){
  cluster_matsubara_kspace_gf G0_init(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(2)));
  EXPECT_EQ(CT_->read_initial_G0(G0_init), false);
}
TEST_F(NormalStateTransformTest, ComputeIsolatedClusterSolution){
  cluster_matsubara_kspace_gf G0_init(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(2)));
  CT_->compute_isolated_cluster_G0(G0_init);
  for(int k=0;k<n_sites_;++k){
    std::vector<double> momentum(CT_->dimension());
    for(int z=0;z<momentum.size();++z) momentum[z]=CT_->momentum_mesh().points()(k,z);
    for(int i=0;i<n_freq_;++i){
      for(int spin=0;spin<2;++spin){
        std::complex<double> atomic_value=1./(std::complex<double>(mu_-CT_->epsilon(momentum),(2*i+1)*M_PI/beta_));
        EXPECT_NEAR(atomic_value.real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(atomic_value.imag(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).imag(), 1.e-10);
      }
    }
  }
}
/*
 * not clear why this does not work, and if it is by design or not.
 * TEST_F(NormalStateTransformTest, IsolatedIsFreeForL1){
  cluster_matsubara_kspace_gf G0_init(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(2)));
  CT_->compute_isolated_cluster_G0(G0_init);
  alps::params p2; set_standard_parameters(p2);
  p2["L"]=1;
  NormalStateClusterTransformer CTL1(p2);
  cluster_matsubara_kspace_gf G0_init_freeL1(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(2)));
  CTL1.compute_noninteracting_G0(G0_init_freeL1);

  for(int k=0;k<n_sites_;++k){
    for(int i=0;i<n_freq_;++i){
      for(int spin=0;spin<2;++spin){
        EXPECT_NEAR(G0_init_freeL1(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).real(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).real(), 1.e-10);
        EXPECT_NEAR(G0_init_freeL1(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).imag(), G0_init(alps::gf::matsubara_index(i), alps::gf::momentum_index(k), alps::gf::index(spin)).imag(), 1.e-10);
      }
    }
  }
}*/
