// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"
#include "cluster_helper.h"
#include "alps/params.hpp"
#include "alps/hdf5.hpp"
#include <cmath>
#include <cstdio>
#include <string>

#define TWODCLUSTERS1 "Betts2D-1A","Betts2D-4A","Betts2D-8A","Betts2D-8C","Betts2D-9A","Betts2D-9B","Betts2D-10A","Betts2D-10B","Betts2D-11A","Betts2D-12A","Betts2D-12B","Betts2D-12D","Betts2D-13A","Betts2D-14A","Betts2D-14B","Betts2D-15A","Betts2D-15B","Betts2D-16A","Betts2D-16B"
//,"Betts2D-16C","Betts2D-16D","Betts2D-16E","Betts2D-17A","Betts2D-17B","Betts2D-18A","Betts2D-18B","Betts2D-18C","Betts2D-18D","Betts2D-19A"
#define TWODCLUSTERS2 "Betts2D-19B","Betts2D-20A","Betts2D-20B","Betts2D-20C","Betts2D-20E"
//,"Betts2D-21A","Betts2D-21B","Betts2D-21C","Betts2D-22A","Betts2D-22B","Betts2D-22D","Betts2D-23A","Betts2D-24A","Betts2D-24B","Betts2D-24C","Betts2D-24D","Betts2D-24E","Betts2D-24F","Betts2D-24G","Betts2D-25A"
#define TWODCLUSTERS3 "Betts2D-25B","Betts2D-26A","Betts2D-26B","Betts2D-26C"
//,"Betts2D-27A","Betts2D-27B","Betts2D-28A","Betts2D-28B","Betts2D-28C","Betts2D-28D","Betts2D-29A","Betts2D-29B","Betts2D-30A","Betts2D-30B","Betts2D-30C","Betts2D-30D","Betts2D-30G","Betts2D-31A","Betts2D-31B","Betts2D-32A","Betts2D-32B","Betts2D-32C","Betts2D-32D","Betts2D-32E"
#define THREEDCLUSTERS1 "Kent3D-8A","Kent3D-8B","Kent3D-8C"

//,"Kent3D-10A","Kent3D-10B","Kent3D-12A","Kent3D-12B","Kent3D-14A","Kent3D-14B","Kent3D-16C","Kent3D-16D","Kent3D-18A","Kent3D-18B","Kent3D-20A","Kent3D-20B","Kent3D-22A","Kent3D-22B","Kent3D-24C","Kent3D-24D","Kent3D-26A","Kent3D-26B","Kent3D-28A","Kent3D-28B","Kent3D-30A","Kent3D-30B","Kent3D-32A","Kent3D-32B","Kent3D-34A","Kent3D-34B","Kent3D-36A","Kent3D-36B","Kent3D-38A","Kent3D-38B","Kent3D-40A","Kent3D-40B","Kent3D-42A","Kent3D-42B","Kent3D-44A","Kent3D-44B","Kent3D-46A"
  void set_other_parameters(alps::params &p){
    p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
    p["MU"] = 0.;
    p["t"] = 1.;
    p["tprime"] = 0.;
    p["BETA"] = 7.;
    p["FLAVORS"] = 2;
    p["NMATSUBARA"] = 1024;
    p["N"] = 1024;
    p["dca.AFM_SELFCONSISTENCY"]=false; //note that we need to know at this point if we're paramagnetic or not. Not sure this is the right place for it.
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

// Setting up the same cluster twice gives the same tables.
TEST_P(ClusterHelperTest, SetupIsReproducible){
  alps::params p2;
  p2["dca.LATTICE"]=GetParam();
  p2["dca.SITES"]=extract_unit_cell_volume_from_name(GetParam());
  p2["dca.L"] = L_;
  set_other_parameters(p2);
  ClusterHelper CH2(p2);

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

// Periodicity: the number of hops along each lattice direction until the same
// cluster site is reached again. The expected values are the smallest n for
// which n e_i is an integer combination of the superlattice vectors. Betts3D-16E,
// -24E and Fuchs3D-50B used to hang: their periodicity was searched for along the
// reciprocal basis vectors, which never ends for a vector without a positive
// component.
struct periodicity_case { const char* lattice; int sites; std::vector<int> periodicity; };
class ClusterPeriodicity : public ::testing::TestWithParam<periodicity_case> {};
TEST_P(ClusterPeriodicity, HopsUntilTheSameSite){
  const periodicity_case& c = GetParam();
  alps::params p;
  p["dca.LATTICE"]=c.lattice;
  p["dca.SITES"]=c.sites;
  p["dca.L"]=4;
  set_other_parameters(p);
  ClusterHelper CH(p);
  for (std::size_t i = 0; i < c.periodicity.size(); ++i)
    EXPECT_EQ(CH.periodicity(i), c.periodicity[i]) << "direction " << i;
}
INSTANTIATE_TEST_SUITE_P(Clusters, ClusterPeriodicity, ::testing::Values(
    periodicity_case{"Betts2D-4A", 4, {2, 2}},
    periodicity_case{"Betts2D-8A", 8, {4, 4}},
    periodicity_case{"Betts2D-8C", 8, {8, 4}},
    periodicity_case{"Betts2D-10A", 10, {10, 5}},
    periodicity_case{"dimer", 2, {2, 1}},
    periodicity_case{"Betts3D-16E", 16, {4, 4, 8}},
    periodicity_case{"Betts3D-24E", 24, {6, 4, 12}},
    periodicity_case{"Fuchs3D-50B", 50, {50, 10, 50}}),
    [](const ::testing::TestParamInfo<periodicity_case>& info) {
      std::string name = info.param.lattice;
      for (char& ch : name) if (!std::isalnum(static_cast<unsigned char>(ch))) ch = '_';
      return name;
    });

// dca.L must be even and at least 2: the patch momenta are integrated with
// Simpson's rule, and for odd L its weights did not add up to L^d, by which
// the coarse-graining divides. Betts2D-4A at L = 9 gave sum(w)/L^2 = 1.075
// and <eps^2> wrong by 0.14 (#441).
TEST(ClusterHelper, RefusesOddOrNonpositiveL){
  for (int L : {-4, 0, 1, 5, 9, 33}) {
    alps::params p;
    p["dca.LATTICE"]="Betts2D-4A";
    p["dca.SITES"]=4;
    p["dca.L"]=L;
    set_other_parameters(p);
    try {
      ClusterHelper CH(p);
      ADD_FAILURE() << "dca.L = " << L << " was accepted";
    } catch (std::invalid_argument const & e) {
      EXPECT_NE(std::string(e.what()).find("dca.L = " + std::to_string(L) + ": "), std::string::npos) << e.what();
      EXPECT_NE(std::string(e.what()).find("must be even and at least 2"), std::string::npos) << e.what();
    }
  }
}

// With an even L the weights of each patch add up to L^d, and the patches
// reproduce the Brillouin-zone average of eps(k) = -2 sum_d cos(k_d):
// <eps> = 0 and <eps^2> = 2d.
struct weight_case { const char* lattice; int sites; int L; };
class ClusterIntegrationWeights : public ::testing::TestWithParam<weight_case> {};
TEST_P(ClusterIntegrationWeights, NormalizedAndExactForTheDispersion){
  const weight_case& c = GetParam();
  alps::params p;
  p["dca.LATTICE"]=c.lattice;
  p["dca.SITES"]=c.sites;
  p["dca.L"]=c.L;
  set_other_parameters(p);
  ClusterHelper CH(p);
  const int d = CH.dimension();
  double wsum = 0, e1 = 0, e2 = 0;
  for (int kc = 0; kc < CH.n_cluster_sites(); ++kc)
    for (int kl = 0; kl < CH.n_lattice_momenta(); ++kl) {
      double e = 0;
      for (int i = 0; i < d; ++i)
        e -= 2 * std::cos(CH.lattice_momentum(kl, i) + CH.cluster_momentum(kc, i));
      const double w = CH.lattice_momentum(kl, d);
      if (kc == 0) wsum += w;
      e1 += w * e;
      e2 += w * e * e;
    }
  const double norm = double(CH.n_cluster_cells()) * CH.n_cluster_sites();
  EXPECT_NEAR(wsum / CH.n_cluster_cells(), 1., 1e-12);
  EXPECT_NEAR(e1 / norm, 0., 1e-12);
  EXPECT_NEAR(e2 / norm, 2. * d, 1e-10);
}
INSTANTIATE_TEST_SUITE_P(Clusters, ClusterIntegrationWeights, ::testing::Values(
    weight_case{"Betts2D-4A", 4, 8},
    weight_case{"Betts2D-8A", 8, 16},
    weight_case{"Betts2D-10A", 10, 2},
    weight_case{"Betts3D-8A", 8, 4}),
    [](const ::testing::TestParamInfo<weight_case>& info) {
      std::string name = std::string(info.param.lattice) + "_L" + std::to_string(info.param.L);
      for (char& ch : name) if (!std::isalnum(static_cast<unsigned char>(ch))) ch = '_';
      return name;
    });
