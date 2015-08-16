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
