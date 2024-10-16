#include "gtest/gtest.h"
#include "cluster_transform.h"
#include "alps/params.hpp"
#include<iostream>

class ClusterTransformTest : public ::testing::Test{
protected:
  virtual void SetUp() {
    alps::params p;
    L_=20;
    p["dca.LATTICE_LIBRARY"] = "../../libcluster/cluster.xml";
    p["dca.LATTICE_DUMP"] = "./cluster.h5";
    p["dca.LATTICE"] = "Betts2D-4A";
    p["dca.SITES"]=4;
    p["dca.L"] = L_;
    p["MU"] = 0.;
    p["t"] = 1.;
    p["U"] = 3.14159;
    p["V"] = 1.2;
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
    CT_=new NormalStateClusterTransformer(p);

    n_freq_=20;
    n_tau_=40;
    n_flavors_=2;
    n_sites_=4;
    beta_=40;
    srand48(0);

  }

  virtual void TearDown() { delete CT_;}
  void initialize_to_random(cluster_matsubara_kspace_gf &G) const{
    for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
      for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
        for (alps::gf::index flavor(0); flavor<G.mesh3().extent(); ++flavor) {
          G(w,k,flavor)=std::complex<double>(drand48(),drand48());
        }
      }
    }
  }
  void initialize_to_random(cluster_itime_kspace_gf &G) const{
    for (alps::gf::itime_index w(0); w<G.mesh1().extent(); ++w) {
      for (alps::gf::momentum_index k(0); k<G.mesh2().extent(); ++k) {
        for (alps::gf::index flavor(0); flavor<G.mesh3().extent(); ++flavor) {
          G(w,k,flavor)=drand48();
        }
      }
    }
  }
  void initialize_to_random(cluster_matsubara_realspace_gf &G) const{
    for (alps::gf::matsubara_index w(0); w<G.mesh1().extent(); ++w) {
      for (alps::gf::real_space_index r1(0); r1<G.mesh2().extent(); ++r1) {
        for (alps::gf::real_space_index r2(0); r2<G.mesh3().extent(); ++r2) {
          for (alps::gf::index flavor(0); flavor<G.mesh4().extent(); ++flavor) {
            G(w,r1,r2,flavor)=std::complex<double>(drand48(),drand48());
          }
        }
      }
    }
  }
  void initialize_to_random(cluster_itime_realspace_gf &G) const{
    for (alps::gf::itime_index w(0); w<G.mesh1().extent(); ++w) {
      for (alps::gf::real_space_index r1(0); r1<G.mesh2().extent(); ++r1) {
        for (alps::gf::real_space_index r2(0); r2<G.mesh3().extent(); ++r2) {
          for (alps::gf::index flavor(0); flavor<G.mesh4().extent(); ++flavor) {
            G(w,r1,r2,flavor)=drand48();
          }
        }
      }
    }
  }
  int L_;
  NormalStateClusterTransformer *CT_;
  int n_freq_;
  int n_flavors_;
  int n_sites_;
  int n_tau_;
  double beta_;
};

TEST_F(ClusterTransformTest, MomentumMeshMakesSense){
  alps::gf::momentum_index_mesh mmesh=CT_->momentum_mesh();
  alps::gf::momentum_index_mesh::container_type four_site_mesh(4, 2);
  four_site_mesh(0,0)=0.; four_site_mesh(0, 1)=0.;
  four_site_mesh(1,0)=M_PI; four_site_mesh(1,1)=M_PI;
  four_site_mesh(2,0)=0.; four_site_mesh(2,1)=M_PI;
  four_site_mesh(3,0)=M_PI; four_site_mesh(3,1)=0.;

  for(int i=0;i<CT_->n_cluster_sites();++i){
    for(int j=0;j<CT_->dimension(); ++j){
      EXPECT_NEAR(four_site_mesh(i, j), mmesh.points()(i, j), 1.e-10);
    }
  }
}
TEST_F(ClusterTransformTest, RealSpaceMeshMakesSense){
  alps::gf::real_space_index_mesh rmesh=CT_->real_space_mesh();
  alps::gf::real_space_index_mesh::container_type four_site_mesh(4, 2);
  four_site_mesh(0, 0)=0.; four_site_mesh(0, 1)=0.;
  four_site_mesh(1, 0)=0.; four_site_mesh(1, 1)=1.;
  four_site_mesh(2, 0)=1.; four_site_mesh(2, 1)=1.;
  four_site_mesh(3, 0)=1.; four_site_mesh(3, 1)=0.;

  for(int i=0;i<CT_->n_cluster_sites();++i){
    for(int j=0;j<CT_->dimension(); ++j){
      EXPECT_NEAR(four_site_mesh(i, j), rmesh.points()(i, j), 1.e-10);
    }
  }
}
TEST_F(ClusterTransformTest, FourierToKAndBack){
  cluster_matsubara_kspace_gf G_ks(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_,n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(n_flavors_)));
  initialize_to_random(G_ks);
  cluster_matsubara_realspace_gf G_rs=CT_->transform_into_real_space(G_ks);
  cluster_matsubara_kspace_gf G_ks2=CT_->transform_into_k_space(G_rs);
  cluster_matsubara_realspace_gf G_rs2=CT_->transform_into_real_space(G_ks2);

  EXPECT_NEAR((G_rs-G_rs2).norm(), 0., 1.e-10);
  EXPECT_NEAR((G_ks-G_ks2).norm(), 0., 1.e-10);
}
TEST_F(ClusterTransformTest, FourierToKAndBackItime){
  cluster_itime_kspace_gf G_ks(alps::gf::itime_k_sigma_gf(alps::gf::itime_mesh(beta_,n_tau_), CT_->momentum_mesh(), alps::gf::index_mesh(n_flavors_)));
  initialize_to_random(G_ks);
  cluster_itime_realspace_gf G_rs=CT_->transform_into_real_space(G_ks);
  cluster_itime_kspace_gf G_ks2=CT_->transform_into_k_space(G_rs);
  cluster_itime_realspace_gf G_rs2=CT_->transform_into_real_space(G_ks2);

  EXPECT_NEAR((G_rs-G_rs2).norm(), 0., 1.e-10);
  EXPECT_NEAR((G_ks-G_ks2).norm(), 0., 1.e-10);
}
TEST_F(ClusterTransformTest, SpinSymmetrySymmetrizes){
  cluster_matsubara_kspace_gf G_w(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_,n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(n_flavors_)));
  cluster_itime_kspace_gf G_t(alps::gf::itime_k_sigma_gf(alps::gf::itime_mesh(beta_,n_tau_), CT_->momentum_mesh(), alps::gf::index_mesh(n_flavors_)));

  initialize_to_random(G_w);
  initialize_to_random(G_t);
  cluster_matsubara_kspace_gf G_w_bup(G_w);
  cluster_itime_kspace_gf G_t_bup(G_t);

  CT_->spin_symmetrize(G_w);
  CT_->spin_symmetrize(G_t);
  for (alps::gf::matsubara_index om=alps::gf::matsubara_index(0); om<G_w.mesh1().extent(); ++om) {
    for (alps::gf::momentum_index k=alps::gf::momentum_index(0); k<G_w.mesh2().extent(); ++k) {
      EXPECT_NEAR(G_w(om,k,alps::gf::index(0)).real(),0.5*(G_w_bup(om,k,alps::gf::index(0))+G_w_bup(om,k,alps::gf::index(1))).real(), 1.e-10);
      EXPECT_NEAR(G_w(om,k,alps::gf::index(0)).imag(),0.5*(G_w_bup(om,k,alps::gf::index(0))+G_w_bup(om,k,alps::gf::index(1))).imag(), 1.e-10);
    }
  }

  for (alps::gf::itime_index om=alps::gf::itime_index(0); om<G_t.mesh1().extent(); ++om) {
    for (alps::gf::momentum_index k=alps::gf::momentum_index(0); k<G_w.mesh2().extent(); ++k) {
      EXPECT_NEAR(G_t(om,k,alps::gf::index(0)),0.5*(G_t_bup(om,k,alps::gf::index(0))+G_t_bup(om,k,alps::gf::index(1))), 1.e-10);
    }
  }

}

TEST_F(ClusterTransformTest, FourSiteKSpaceSymmetry){
  cluster_matsubara_kspace_gf G_w(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta_,n_freq_), CT_->momentum_mesh(), alps::gf::index_mesh(n_flavors_)));
  cluster_itime_kspace_gf G_t(alps::gf::itime_k_sigma_gf(alps::gf::itime_mesh(beta_,n_tau_), CT_->momentum_mesh(), alps::gf::index_mesh(n_flavors_)));

  initialize_to_random(G_w);
  initialize_to_random(G_t);
  cluster_matsubara_kspace_gf G_w_bup(G_w);
  cluster_itime_kspace_gf G_t_bup(G_t);

  CT_->symmetrize_k_space(G_w);
  CT_->symmetrize_k_space(G_t);

  for (alps::gf::matsubara_index om=alps::gf::matsubara_index(0); om<G_w.mesh1().extent(); ++om) {
    for (alps::gf::index f=alps::gf::index(0); f<G_w.mesh3().extent(); ++f) {
      EXPECT_NEAR(0,std::abs(G_w(om,alps::gf::momentum_index(3),f)-G_w(om,alps::gf::momentum_index(2),f)), 1.e-10);
    }
  }
  for (alps::gf::itime_index om=alps::gf::itime_index(0); om<G_t.mesh1().extent(); ++om) {
    for (alps::gf::index f=alps::gf::index(0); f<G_t.mesh3().extent(); ++f) {
      EXPECT_NEAR(0,std::abs(G_t(om,alps::gf::momentum_index(3),f)-G_t(om,alps::gf::momentum_index(2),f)), 1.e-10);
    }
  }
}
TEST_F(ClusterTransformTest, FourSiteRealSpaceSymmetry){
  cluster_matsubara_realspace_gf G_w(alps::gf::omega_r1_r2_sigma_gf(alps::gf::matsubara_positive_mesh(beta_,n_freq_), CT_->real_space_mesh(),CT_->real_space_mesh(), alps::gf::index_mesh(n_flavors_)));
  cluster_itime_realspace_gf G_t(alps::gf::itime_r1_r2_sigma_gf(alps::gf::itime_mesh(beta_,n_tau_), CT_->real_space_mesh(),CT_->real_space_mesh(), alps::gf::index_mesh(n_flavors_)));

  initialize_to_random(G_w);
  initialize_to_random(G_t);
  cluster_matsubara_realspace_gf G_w_bup(G_w);
  cluster_itime_realspace_gf G_t_bup(G_t);

  CT_->symmetrize_real_space(G_w);
  CT_->symmetrize_real_space(G_t);
  for (alps::gf::matsubara_index om=alps::gf::matsubara_index(0); om<G_w.mesh1().extent(); ++om) {
    for (alps::gf::index f=alps::gf::index(0); f<G_w.mesh4().extent(); ++f) {
      EXPECT_NEAR(0,std::abs(G_w(om,alps::gf::real_space_index(0),alps::gf::real_space_index(1),f)-G_w(om,alps::gf::real_space_index(1),alps::gf::real_space_index(0),f)), 1.e-10);
    }
  }
  for (alps::gf::itime_index om=alps::gf::itime_index(0); om<G_t.mesh1().extent(); ++om) {
    for (alps::gf::index f=alps::gf::index(0); f<G_t.mesh4().extent(); ++f) {
      EXPECT_NEAR(0,std::abs(G_t(om,alps::gf::real_space_index(0),alps::gf::real_space_index(1),f)-G_t(om,alps::gf::real_space_index(1),alps::gf::real_space_index(0),f)), 1.e-10);
    }
  }
}
