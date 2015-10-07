#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"

class ClusterTransformTest : public ::testing::Test{
 protected:
  virtual void SetUp() {
    alps::params p;
    L_=20;
    p["LATTICE_LIBRARY"] = "./cluster.xml";
    p["LATTICE_DUMP"] = "./cluster.h5";
    p["LATTICE"] = "Betts2D-4A";
    p["SITES"]=4;
    p["L"] = L_;
    p["AFM_SELFCONSISTENCY"]=false; //note that we need to know at this point if we're paramagnetic or not. Not sure this is the right place for it.
    CT_=new NormalStateClusterTransformer(p);
  }

  virtual void TearDown() { delete CT_;}
  int L_;
  NormalStateClusterTransformer *CT_;
};

TEST_F(ClusterTransformTest, MeshMakesSense){
  alps::gf::momentum_index_mesh mmesh=CT_->momentum_mesh();
  alps::gf::momentum_index_mesh::container_type four_site_mesh(boost::extents[4][2]);
  four_site_mesh[0][0]=0.; four_site_mesh[0][1]=0.;
  four_site_mesh[1][0]=M_PI; four_site_mesh[1][1]=M_PI;
  four_site_mesh[2][0]=0.; four_site_mesh[2][1]=M_PI;
  four_site_mesh[3][0]=M_PI; four_site_mesh[3][1]=0.;

  for(int i=0;i<CT_->n_cluster_sites();++i){
    for(int j=0;j<CT_->dimension(); ++j){
      EXPECT_NEAR(four_site_mesh[i][j], mmesh.points()[i][j], 1.e-10);
    }
  }
}

