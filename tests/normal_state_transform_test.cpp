#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"
#include<iostream>

class NormalStateTransformTest : public ::testing::Test{
protected:
  virtual void SetUp() {
    alps::params p;
    set_standard_parameters(p);
    CT_=new DoubleUnitCellClusterTransformer(p);

    n_freq_=20;
    n_tau_=40;
    n_flavors_=2;
    n_sites_=4;
    beta_=40;
    srand48(0);

  }
  virtual void TearDown() {
    delete CT_;
  }
  void set_standard_parameters(alps::params &p){
    p["dca.L"] = 20;
    p["dca.LATTICE"] = "Betts2D-8A";
    p["dca.SITES"] = 8;
    p["dca.LATTICE_LIBRARY"] = "../libcluster/cluster.xml";
    p["dca.LATTICE_DUMP"] = "./cluster.h5";
    p["MU"] = 0.;
    p["t"] = 1.;
    p["tprime"] = 0.;
    p["t_flavor"] = 0.;
    p["BETA"] = 7.;
    p["FLAVORS"] = 2;
    p["NMATSUBARA"] = 1024;
    p["N"] = 1024;
    p["dca.AFM_SELFCONSISTENCY"]=false; //note that we need to know at this point if we're paramagnetic or not. Not sure this is the right place for it.
    p["dca.PHASE"]=false;
    p["dca.IMBALANCE"]=false;
    p["dca.G0_OMEGA_INPUT"]="/tmp/none";
    p["dca.G0_OMEGA_OUTPUT"]="/tmp/none";
    p["dca.INITIAL_G0"]="FREE";
  }
  int L_;
  DoubleUnitCellClusterTransformer *CT_;
  int n_freq_;
  int n_flavors_;
  int n_sites_;
  int n_tau_;
  double beta_;
};
TEST_F(NormalStateTransformTest, MakeANormalStateTransformer){
 ; //just run Setup/TearDown
}
TEST_F(NormalStateTransformTest, ReadNonexistentInputFile){
  cluster_matsubara_kspace_gf G0_init(omega_k_sigma_gf_with_tail(alps::gf::matsubara_positive_mesh(beta_, n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(2)));
 EXPECT_EQ(CT_->read_initial_G0(G0_init), false);
}
