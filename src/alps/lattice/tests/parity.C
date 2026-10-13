/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/lattice.h>
#include <array>
#include <ostream>
#include <vector>

namespace {
enum class Pattern { checkerboard, undefined, black, stripes };
struct ParityScenario {
  const char* name;
  bool diagonals;
  const char* backbone;
  Pattern pattern;
};
void PrintTo(const ParityScenario& scenario, std::ostream* out) { *out << scenario.name; }
class LatticeParity : public ::testing::TestWithParam<ParityScenario> {};

void expect_periodic_square(const alps::graph_helper<>& lattice, bool diagonals) {
  ASSERT_EQ(lattice.dimension(), 2u);
  ASSERT_EQ(lattice.num_sites(), 16u);
  const int bonds_per_site = diagonals ? 4 : 2;
  ASSERT_EQ(lattice.num_bonds(), 16u * bonds_per_site);
  const auto& graph = lattice.graph();
  for (int site = 0; site < 16; ++site) {
    SCOPED_TRACE(::testing::Message() << "site " << site);
    EXPECT_EQ(lattice.site_type(site), 0);
    EXPECT_EQ(lattice.coordinate(site), (std::vector<double>{double(site / 4), double(site % 4)}));
  }
  // The former XML recorded every endpoint, bond type and oriented displacement.
  // These four unit-cell directions independently specify the same periodic graph.
  const std::array<std::array<int, 2>, 4> directions{{{{0, 1}}, {{1, 0}}, {{1, 1}}, {{1, -1}}}};
  const auto edges = boost::edges(graph);
  int index = 0;
  for (auto edge = edges.first; edge != edges.second; ++edge, ++index) {
    SCOPED_TRACE(::testing::Message() << "bond " << index);
    const int site = index / bonds_per_site;
    const int direction = index % bonds_per_site;
    const auto delta = directions[direction];
    const int x = (site / 4 + delta[0]) % 4;
    const int y = (site % 4 + delta[1] + 4) % 4;
    EXPECT_EQ(get(boost::edge_index, graph, *edge), index);
    EXPECT_EQ(boost::source(*edge, graph), site);
    EXPECT_EQ(boost::target(*edge, graph), x * 4 + y);
    EXPECT_EQ(lattice.bond_type(*edge), direction < 2 ? 0 : 1);
    EXPECT_EQ(lattice.bond_vector(*edge), (std::vector<double>{double(delta[0]), double(delta[1])}));
  }
}

TEST_P(LatticeParity, AppliesBackboneWithoutChangingGraph) {
  const auto& scenario = GetParam();
  alps::Parameters parameters;
  parameters["LATTICE_LIBRARY"] = "lattices.xml";
  parameters["LATTICE"] = scenario.diagonals ? "frustrated square lattice" : "square lattice";
  parameters["L"] = 4;
  if (scenario.backbone) parameters["BACKBONE_TYPES"] = scenario.backbone;
  alps::graph_helper<> lattice(parameters);
  ASSERT_NO_FATAL_FAILURE(expect_periodic_square(lattice, scenario.diagonals));
  EXPECT_EQ(lattice.is_bipartite(), scenario.pattern != Pattern::undefined);

  using Traits = alps::parity_traits<alps::parity_t, alps::coordinate_graph_type>;
  for (int site = 0; site < 16; ++site) {
    SCOPED_TRACE(::testing::Message() << "site " << site);
    auto expected = Traits::black;
    switch (scenario.pattern) {
      case Pattern::checkerboard:
        expected = (site / 4 + site % 4) % 2 ? Traits::white : Traits::black;
        break;
      case Pattern::undefined: expected = Traits::undefined; break;
      case Pattern::black: break;
      case Pattern::stripes: expected = (site / 4) % 2 ? Traits::white : Traits::black; break;
    }
    EXPECT_EQ(get(alps::parity_t(), lattice.graph(), site), expected);
    EXPECT_EQ(lattice.parity(site), expected == Traits::undefined ? 0.0 : expected == Traits::black ? -1.0 : 1.0);
  }
}

// All six original parity.input scenarios; empty and absent backbones differ.
INSTANTIATE_TEST_SUITE_P(HistoricalBackbones, LatticeParity, ::testing::Values(
  ParityScenario{"SquareDefault", false, nullptr, Pattern::checkerboard},
  ParityScenario{"FrustratedDefault", true, nullptr, Pattern::undefined},
  ParityScenario{"EmptyBackbone", true, "", Pattern::black},
  ParityScenario{"NearestNeighborBackbone", true, "0", Pattern::checkerboard},
  ParityScenario{"DiagonalBackbone", true, "1", Pattern::stripes},
  ParityScenario{"AllBondsBackbone", true, "0 1", Pattern::undefined}),
  [](const ::testing::TestParamInfo<ParityScenario>& info) { return info.param.name; });
} // namespace
