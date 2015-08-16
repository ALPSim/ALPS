#include "gtest/gtest.h"
#include "alpsgraph.h"

TEST(ALPSGraphTest, ALPSGraphCreation){
  alps::params p;
  p["LATTICE_LIBRARY"] = "./cluster.xml";
  p["LATTICE"] = "Betts2D-4A";
  p["L"] = "40";
  alps_graph G(p);

  EXPECT_EQ(40*40, G.volume());
  EXPECT_EQ(2, G.dimension());
}
TEST(ALPSGraphTest, BasisVectors){
  alps::params p;
  p["LATTICE_LIBRARY"] = "./cluster.xml";
  p["LATTICE"] = "Betts2D-4A";
  p["L"] = "40";
  alps_graph G(p);
  
  alps_graph::basis_vector_iterator_pair BVIP=G.basis_vectors();
  
  alps_graph::basis_vector basis(G.dimension(), std::vector<double>(2, 0.));
  basis[0][0]=2.;   basis[1][1]=2.;
  
  std::vector<std::vector<double> >::const_iterator my_it=basis.begin();
  for(std::vector<std::vector<double> >::const_iterator it=BVIP.first; it!=BVIP.second;++it){
    for(int i=0;i<G.dimension();++i){
      EXPECT_EQ((*my_it)[i], (*it)[i]);
    }
    ++my_it;
  }
}
TEST(ALPSGraphTest, ReciprocalVectors){
  alps::params p;
  p["LATTICE_LIBRARY"] = "./cluster.xml";
  p["LATTICE"] = "Betts2D-4A";
  p["L"] = "40";
  int size=4;
  alps_graph G(p);
  
  alps_graph::basis_vector_iterator_pair RVIP=G.reciprocal_basis_vectors();
  
  alps_graph::basis_vector basis           (G.dimension(), std::vector<double>(2, 0.));
  alps_graph::basis_vector reciprocal_basis(G.dimension(), std::vector<double>(2, 0.));
  basis[0][0]=2.;   basis[1][1]=2.;
  
  for(int i=0;i<G.dimension();++i){
    for(int j=0;j<G.dimension();++j){
      reciprocal_basis[i][j]=basis[i][j]*2*M_PI/size;
    }
  }

  alps_graph::basis_vector_iterator my_it=reciprocal_basis.begin();
  for(alps_graph::basis_vector_iterator it=RVIP.first; it!=RVIP.second;++it){
    for(int i=0;i<G.dimension();++i){
      EXPECT_EQ((*my_it)[i], (*it)[i]);
    }
    ++my_it;
  }
}

#define TWODCLUSTERS1 "Betts2D-1A","Betts2D-4A","Betts2D-8A","Betts2D-8C","Betts2D-9A","Betts2D-9B","Betts2D-10A","Betts2D-10B","Betts2D-11A","Betts2D-12A","Betts2D-12B","Betts2D-12D","Betts2D-13A","Betts2D-14A","Betts2D-14B","Betts2D-15A","Betts2D-15B","Betts2D-16A","Betts2D-16B","Betts2D-16C","Betts2D-16D","Betts2D-16E","Betts2D-17A","Betts2D-17B","Betts2D-18A","Betts2D-18B","Betts2D-18C","Betts2D-18D","Betts2D-19A"
#define TWODCLUSTERS2 "Betts2D-19B","Betts2D-20A","Betts2D-20B","Betts2D-20C","Betts2D-20E","Betts2D-21A","Betts2D-21B","Betts2D-21C","Betts2D-22A","Betts2D-22B","Betts2D-22D","Betts2D-23A","Betts2D-24A","Betts2D-24B","Betts2D-24C","Betts2D-24D","Betts2D-24E","Betts2D-24F","Betts2D-24G","Betts2D-25A"
#define TWODCLUSTERS3 "Betts2D-25B","Betts2D-26A","Betts2D-26B","Betts2D-26C","Betts2D-27A","Betts2D-27B","Betts2D-28A","Betts2D-28B","Betts2D-28C","Betts2D-28D","Betts2D-29A","Betts2D-29B","Betts2D-30A","Betts2D-30B","Betts2D-30C","Betts2D-30D","Betts2D-30G","Betts2D-31A","Betts2D-31B","Betts2D-32A","Betts2D-32B","Betts2D-32C","Betts2D-32D","Betts2D-32E","Betts2D-98A"
#define THREEDCLUSTERS1 "Kent3D-8A","Kent3D-8B","Kent3D-8C","Kent3D-10A","Kent3D-10B","Kent3D-12A","Kent3D-12B","Kent3D-14A","Kent3D-14B","Kent3D-16C","Kent3D-16D","Kent3D-18A","Kent3D-18B","Kent3D-20A","Kent3D-20B","Kent3D-22A","Kent3D-22B","Kent3D-24C","Kent3D-24D","Kent3D-26A","Kent3D-26B","Kent3D-28C","Kent3D-28D","Kent3D-30A","Kent3D-30B","Kent3D-32C","Kent3D-32D","Kent3D-34A","Kent3D-34B","Kent3D-36A","Kent3D-36B","Kent3D-38A","Kent3D-38B","Kent3D-40A","Kent3D-40B","Kent3D-42A","Kent3D-42B","Kent3D-44A","Kent3D-44B","Kent3D-46A"

#define THREEDCLUSTERS2 "Kent3D-46B","Kent3D-48A","Kent3D-48B","Kent3D-50A","Kent3D-50B","Kent3D-52A","Kent3D-52B","Kent3D-54A","Kent3D-54B","Kent3D-56A","Kent3D-56B","Kent3D-58A","Kent3D-58B","Kent3D-60A","Kent3D-60B","Kent3D-62A","Kent3D-62B","Kent3D-64A","Kent3D-64B","Kent3D-66A","Kent3D-66B","Kent3D-68A","Kent3D-68B","Kent3D-70A","Kent3D-70B","Kent3D-72A","Kent3D-72B","Kent3D-74A","Kent3D-74B","Kent3D-76A","Kent3D-76B","Kent3D-80A","Kent3D-80B","Kent3D-84A","Kent3D-84B","Kent3D-88A","Kent3D-88B","Kent3D-90A","Kent3D-90B","Kent3D-96A"
#define THREEDCLUSTERS3 "Kent3D-100A","Kent3D-102A","Kent2D-4A","Kent2D-8A","Kent2D-8B","Kent2D-10A","Kent2D-10B","Kent2D-12A","Kent2D-12B","Kent2D-14A","Kent2D-14B","Kent2D-16A","Kent2D-16B","Kent2D-18A","Kent2D-18B","Kent2D-20A","Kent2D-22A","Kent2D-22B","Kent2D-24A","Kent2D-24B","Kent2D-24C","Kent2D-26A","Kent2D-26B","Kent2D-28A","Kent2D-28B","Kent2D-30A","Kent2D-30B","Kent2D-32A","Kent2D-32B","Kent2D-34A","Kent2D-34B","Kent2D-36A","Kent2D-36B","Kent2D-38A","Kent2D-38B","Kent2D-38C","Kent2D-40A","Kent2D-40B","Kent2D-40C","Kent2D-42A","Kent2D-42B","Kent2D-44A"
#define THREEDCLUSTERS4 "Kent2D-44B","Kent2D-46A","Kent2D-46B","Kent2D-48A","Kent2D-48B","Kent2D-50A","Kent2D-50B","Kent2D-52A","Kent2D-52B","Kent2D-54A","Kent2D-54B","Kent2D-56A","Kent2D-56B","Kent2D-58A","Kent2D-58B","Kent2D-60A","Kent2D-60B","Kent2D-62A","Kent2D-62B","Kent2D-64A","Kent2D-64B","Kent2D-66A","Kent2D-66B","Kent2D-68A","Kent2D-68B","Kent2D-70A","Kent2D-70B","Kent2D-72A","Kent2D-72B","Kent2D-74A","Kent2D-74B","Kent2D-76A","Kent2D-76B","Kent2D-78A","Kent2D-78B","Kent2D-80A","Kent2D-80B","Kent2D-180A","Kent2D-182A","Kent2D-182B","Kent2D-184A","Kent2D-184B"
#define THREEDCLUSTERS5 "Kent2D-186A","Kent2D-186B","Kent2D-188A","Kent2D-188B","Kent2D-190A","Kent2D-190B","Kent2D-192A","Kent2D-192B","Kent2D-194A","Kent2D-194B","Kent2D-196A","Kent2D-196B","Kent2D-198A","Kent2D-198B","Kent2D-200A","Kent2D-200B","Kent2D-202A","Kent2D-202B","Kent2D-204A","Kent2D-204B","Kent2D-206A","Kent2D-206B","Kent2D-208A","Kent2D-208B","Kent2D-210A","Kent2D-210B","Kent2D-212A","Kent2D-212B","Kent2D-214A","Kent2D-214B","Kent2D-216A","Kent2D-216B","Kent2D-218A","Kent2D-218B","Kent2D-220A","Kent2D-220B","Kent2D-222A","Kent2D-222B","Kent2D-224A"
#define THREEDCLUSTERS6 "Kent2D-224B","Kent2D-226A","Kent2D-226B","Kent2D-228A","Kent2D-228B","Kent2D-230A","Kent2D-230B","Kent2D-232A","Kent2D-232B","Kent2D-234A","Kent2D-234B","Kent2D-236A","Kent2D-236B","Kent2D-238A","Kent2D-238B","Kent2D-240A","Kent2D-240B","Kent2D-242A","Kent2D-242B","Kent2D-244A","Kent2D-244B","Kent2D-246A","Kent2D-246B","Kent2D-248A","Kent2D-248B","Kent2D-250A","Kent2D-250B","Kent2D-252A","Kent2D-252B","Kent2D-254A","Kent2D-254B","Kent2D-256A","Kent2D-256B","Kent2D-258A","Kent2D-258B","Kent2D-260A","Kent2D-260B","Kent2D-262A","Kent2D-262B"
#define THREEDCLUSTERS7 "Kent2D-264A","Kent2D-264B","Kent2D-266A","Kent2D-266B","Kent2D-268A","Kent2D-268B","Kent2D-270A","Kent2D-270B","Kent2D-272A","Kent2D-272B","Kent2D-274A","Kent2D-274B","Kent2D-276A","Kent2D-276B","Kent2D-278A","Kent2D-278B","Kent2D-280A","Kent2D-280B","Kent2D-282A","Kent2D-282B","Kent2D-284A","Kent2D-284B","Kent2D-286A","Kent2D-286B","Kent2D-288A","Kent2D-288B","Kent2D-290A","Kent2D-290B","Kent2D-292A","Kent2D-292B","Kent2D-294A","Kent2D-294B","Kent2D-296A","Kent2D-296B","Kent2D-298A","Kent2D-298B","Kent2D-300A","Kent2D-300B","Kent2D-302A"
#define THREEDCLUSTERS8 "Kent2D-302B","Kent2D-304A","Kent2D-304B","Kent2D-306A","Kent2D-306B","Kent2D-308A","Kent2D-308B","Kent2D-310A","Kent2D-310B","Kent2D-312A","Kent2D-312B","Kent2D-314A","Kent2D-314B","Kent2D-316A","Kent2D-316B","Kent2D-318A","Kent2D-318B","Kent2D-320A","Kent2D-320B","Kent2D-322A","Kent2D-322B","Kent2D-324A","Kent2D-324B","Kent2D-326A","Kent2D-326B","Kent2D-328A","Kent2D-328B","Kent2D-330A","Kent2D-330B","Kent2D-332A","Kent2D-332B","Kent2D-334A","Kent2D-334B","Kent2D-336A","Kent2D-336B","Kent2D-338A","Kent2D-338B","Kent2D-340A","Kent2D-340B"
#define THREEDCLUSTERS9 "Kent2D-342A","Kent2D-342B","Kent2D-344A","Kent2D-344B","Kent2D-346A","Kent2D-346B","Kent2D-348A","Kent2D-348B","Kent2D-350A","Kent2D-350B","Kent2D-352A","Kent2D-352B","Kent2D-354A","Kent2D-354B","Kent2D-356A","Kent2D-356B","Kent2D-358A","Kent2D-358B","Kent2D-360A","Kent2D-360B","Kent2D-362A","Kent2D-362B","Kent2D-364A","Kent2D-364B","Kent2D-366A","Kent2D-366B","Kent2D-368A","Kent2D-368B","Kent2D-370A","Kent2D-370B","Kent2D-372A","Kent2D-372B","Kent2D-374A","Kent2D-376A","Kent2D-378A","Kent2D-384A","Kent2D-386A","Kent2D-386B"
#define THREEDCLUSTERS10 "Kent2D-388A","Kent2D-388B","Kent2D-390A","Kent2D-390B","Kent2D-392A","Kent2D-392B","Kent2D-394A","Kent2D-394B","Kent2D-396A","Kent2D-396B","Kent2D-400A","Kent2D-400B","Kent2D-402A","Kent2D-416A","Kent2D-418A","Kent2D-418B","Kent2D-420A","Kent2D-420B","Kent2D-422A","Kent2D-422B","Kent2D-432A","Kent2D-448A","Kent2D-450A","Kent2D-450B","Kent2D-464A","Kent2D-480A","Kent2D-496A","Kent2D-512A"

class ClusterTest : public ::testing::TestWithParam<const char*> {
 protected:
  virtual void SetUp() {
    alps::params p;
    L_=2;
    p["LATTICE_LIBRARY"] = "./cluster.xml";
    p["LATTICE"] = GetParam();
    p["L"] = L_;
    G_=new alps_graph(p);
    extract_unit_cell_volume_from_name(GetParam());
  }
  void extract_unit_cell_volume_from_name(const char* name){
    std::string name_str(name);
    int dashpos=name_str.find('-');
    if(dashpos==std::string::npos) throw std::invalid_argument("could not find location of - in name.");
    name_str.erase(0, dashpos+1);
    std::stringstream n(name_str);
    n>>unit_cell_volume_extracted_;
  }

  virtual void TearDown() { delete G_;}
  int unit_cell_volume_extracted_;
  int L_;
  alps_graph *G_;
};

TEST_P(ClusterTest, TwoDClusterCreation){
  EXPECT_EQ(std::pow((double)L_,G_->dimension()), G_->volume());
  EXPECT_EQ(unit_cell_volume_extracted_, G_->unit_cell_volume());
}
INSTANTIATE_TEST_CASE_P(ClusterCreation2D1 , ClusterTest,::testing::Values(TWODCLUSTERS1));
INSTANTIATE_TEST_CASE_P(ClusterCreation2D2 , ClusterTest,::testing::Values(TWODCLUSTERS2));
INSTANTIATE_TEST_CASE_P(ClusterCreation2D3 , ClusterTest,::testing::Values(TWODCLUSTERS3));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D1 , ClusterTest,::testing::Values(THREEDCLUSTERS1));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D2 , ClusterTest,::testing::Values(THREEDCLUSTERS2));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D3 , ClusterTest,::testing::Values(THREEDCLUSTERS3));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D4 , ClusterTest,::testing::Values(THREEDCLUSTERS4));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D5 , ClusterTest,::testing::Values(THREEDCLUSTERS5));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D6 , ClusterTest,::testing::Values(THREEDCLUSTERS6));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D7 , ClusterTest,::testing::Values(THREEDCLUSTERS7));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D8 , ClusterTest,::testing::Values(THREEDCLUSTERS8));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D9 , ClusterTest,::testing::Values(THREEDCLUSTERS9));
INSTANTIATE_TEST_CASE_P(ClusterCreation3D10, ClusterTest,::testing::Values(THREEDCLUSTERS10));