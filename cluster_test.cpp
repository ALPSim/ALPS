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
