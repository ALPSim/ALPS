// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"
#include "cluster_helper.h"
#include "alps/params.hpp"

#define TWODCLUSTERS1 "Betts2D-1A","Betts2D-4A","Betts2D-8A","Betts2D-8C","Betts2D-9A","Betts2D-9B","Betts2D-10A","Betts2D-10B","Betts2D-11A","Betts2D-12A","Betts2D-12B","Betts2D-12D","Betts2D-13A","Betts2D-14A","Betts2D-14B","Betts2D-15A","Betts2D-15B","Betts2D-16A","Betts2D-16B"
//,"Betts2D-16C","Betts2D-16D","Betts2D-16E","Betts2D-17A","Betts2D-17B","Betts2D-18A","Betts2D-18B","Betts2D-18C","Betts2D-18D","Betts2D-19A"
#define TWODCLUSTERS2 "Betts2D-19B","Betts2D-20A","Betts2D-20B","Betts2D-20C","Betts2D-20E"
//,"Betts2D-21A","Betts2D-21B","Betts2D-21C","Betts2D-22A","Betts2D-22B","Betts2D-22D","Betts2D-23A","Betts2D-24A","Betts2D-24B","Betts2D-24C","Betts2D-24D","Betts2D-24E","Betts2D-24F","Betts2D-24G","Betts2D-25A"
#define TWODCLUSTERS3 "Betts2D-25B","Betts2D-26A","Betts2D-26B","Betts2D-26C"
//,"Betts2D-27A","Betts2D-27B","Betts2D-28A","Betts2D-28B","Betts2D-28C","Betts2D-28D","Betts2D-29A","Betts2D-29B","Betts2D-30A","Betts2D-30B","Betts2D-30C","Betts2D-30D","Betts2D-30G","Betts2D-31A","Betts2D-31B","Betts2D-32A","Betts2D-32B","Betts2D-32C","Betts2D-32D","Betts2D-32E"
#define THREEDCLUSTERS1 "Kent3D-8A","Kent3D-8B","Kent3D-8C"

//,"Kent3D-10A","Kent3D-10B","Kent3D-12A","Kent3D-12B","Kent3D-14A","Kent3D-14B","Kent3D-16C","Kent3D-16D","Kent3D-18A","Kent3D-18B","Kent3D-20A","Kent3D-20B","Kent3D-22A","Kent3D-22B","Kent3D-24C","Kent3D-24D","Kent3D-26A","Kent3D-26B","Kent3D-28C","Kent3D-28D","Kent3D-30A","Kent3D-30B","Kent3D-32C","Kent3D-32D","Kent3D-34A","Kent3D-34B","Kent3D-36A","Kent3D-36B","Kent3D-38A","Kent3D-38B","Kent3D-40A","Kent3D-40B","Kent3D-42A","Kent3D-42B","Kent3D-44A","Kent3D-44B","Kent3D-46A"
  void set_other_parameters(alps::params &p){
    p["dca.LATTICE_LIBRARY"] = "../../libcluster/cluster.xml";
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

class ClusterHelperTest : public ::testing::TestWithParam<const char*> {
 protected:
  virtual void SetUp() {
    alps::params p;
    L_=20;
    p["dca.LATTICE"] = GetParam();
    p["dca.L"] = L_;
    p["dca.SITES"]=extract_unit_cell_volume_from_name(GetParam());
    p["dca.L"] = L_;
    set_other_parameters(p);
    CH_=new ClusterHelper(p);
  }
  int extract_unit_cell_volume_from_name(const char* name){
    int unit_cell_volume_extracted;
    std::string name_str(name);
    int dashpos=name_str.find('-');
    if(dashpos==std::string::npos) throw std::invalid_argument("could not find location of - in name.");
    name_str.erase(0, dashpos+1);
    std::stringstream n(name_str);
    n>>unit_cell_volume_extracted;
    return unit_cell_volume_extracted;
  }

  virtual void TearDown() { delete CH_;}
  int L_;
  ClusterHelper *CH_;
};

TEST_P(ClusterHelperTest, TwoDClusterReadWrite){
  alps::params p2;
  p2["dca.LATTICE"]=GetParam();
  p2["TMPNAME"]="cluster_dump.h5";
  set_other_parameters(p2);
  CH_->write();
  
  p2["dca.SITES"]=extract_unit_cell_volume_from_name(GetParam());
  p2["dca.L"] = L_;
  ClusterHelper CH2(p2);
  CH2.read();
  
  EXPECT_EQ(*CH_, CH2);
}
TEST_P(ClusterHelperTest, TwoDClusterVsStoredALPS){
  alps::params p2;
  p2["dca.LATTICE"]=GetParam();
  p2["dca.TMPNAME"]="cluster_dump.bup.h5";
  p2["dca.SITES"]=extract_unit_cell_volume_from_name(GetParam());
  p2["dca.L"] = L_;
  set_other_parameters(p2);
  ClusterHelper CH2(p2);
  CH2.read();
  //CH2.read(p2, GetParam());

  EXPECT_EQ(*CH_, CH2);
}

INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest2D1 , ClusterHelperTest,::testing::Values(TWODCLUSTERS1));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest2D2 , ClusterHelperTest,::testing::Values(TWODCLUSTERS2));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest2D3 , ClusterHelperTest,::testing::Values(TWODCLUSTERS3));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D1 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS1));
/*INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D2 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS2));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D3 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS3));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D4 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS4));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D5 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS5));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D6 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS6));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D7 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS7));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D8 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS8));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D9 , ClusterHelperTest,::testing::Values(THREEDCLUSTERS9));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D10, ClusterHelperTest,::testing::Values(THREEDCLUSTERS10));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D11, ClusterHelperTest,::testing::Values(THREEDCLUSTERS11));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D12, ClusterHelperTest,::testing::Values(THREEDCLUSTERS12));
INSTANTIATE_TEST_CASE_P(ClusterNoALPSTest3D13, ClusterHelperTest,::testing::Values(THREEDCLUSTERS13));*/

TEST(ClusterHelper, Periodicity8Site){
  alps::params p;
  p["dca.LATTICE"]="Betts2D-8A";
  p["dca.SITES"]=8;
  p["dca.L"] = 20;
  set_other_parameters(p);


  ClusterHelper CH(p);

  EXPECT_EQ(CH.periodicity(0), 4);
  EXPECT_EQ(CH.periodicity(1), 4);
}
TEST(ClusterHelper, Periodicity4Site){
  alps::params p;
  p["dca.LATTICE"]="Betts2D-4A";
  p["dca.SITES"]=4;
  p["dca.L"] = 20;
  p["dca.AFM_SELFCONSISTENCY"]=false;
  set_other_parameters(p);
  ClusterHelper CH(p);

  EXPECT_EQ(CH.periodicity(0), 2);
  EXPECT_EQ(CH.periodicity(1), 2);
}
