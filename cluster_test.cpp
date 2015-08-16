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
class TwoDClusterTest : public ::testing::TestWithParam<const char*> {
 protected:
  virtual void SetUp() {
    alps::params p;
    p["LATTICE_LIBRARY"] = "./cluster.xml";
    p["LATTICE"] = GetParam();
    p["L"] = "40";
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
  alps_graph *G_;
};

TEST_P(TwoDClusterTest, ClusterCreation){
  EXPECT_EQ(40*40, G_->volume());
  EXPECT_EQ(2, G_->dimension());
  EXPECT_EQ(unit_cell_volume_extracted_, G_->unit_cell_volume());
}
INSTANTIATE_TEST_CASE_P(ClusterCreation1, TwoDClusterTest,::testing::Values(TWODCLUSTERS1));
INSTANTIATE_TEST_CASE_P(ClusterCreation2, TwoDClusterTest,::testing::Values(TWODCLUSTERS2));
INSTANTIATE_TEST_CASE_P(ClusterCreation3, TwoDClusterTest,::testing::Values(TWODCLUSTERS3));