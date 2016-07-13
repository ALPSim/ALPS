#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"
#include<iostream>

class DoubleUnitCellTransformTest : public ::testing::Test{
protected:
  virtual void SetUp() {
    alps::params p;
    L_=20;

    CT_=new DoubleUnitCellClusterTransformer(p);

    n_freq_=20;
    n_tau_=40;
    n_flavors_=2;
    n_sites_=4;
    beta_=40;
    srand48(0);

  }
  void set_other_parameters(alps::params &p){
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
TEST_F(DoubleUnitCellTransformTest, MakeADoubleUnitCellTransformer){

}
