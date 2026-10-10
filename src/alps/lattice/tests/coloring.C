/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/lattice.h>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <ostream>
#include <vector>

namespace {
struct ColoringScenario {
  const char* name;
  const char* lattice;
  int length;
  std::size_t count;
  std::vector<std::size_t> colors;
};
void PrintTo(const ColoringScenario& scenario, std::ostream* out) { *out << scenario.name; }
class LatticeColoring : public ::testing::TestWithParam<ColoringScenario> {};

TEST_P(LatticeColoring, PreservesColorsAndSeparatesBondEndpoints) {
  const auto& scenario = GetParam();
  alps::Parameters parameters;
  parameters["LATTICE_LIBRARY"] = "lattices.xml";
  parameters["LATTICE"] = scenario.lattice;
  parameters["L"] = scenario.length;
  alps::graph_helper<> lattice(parameters);
  ASSERT_EQ(lattice.num_sites(), scenario.colors.size());
  std::vector<std::size_t> colors(lattice.num_sites());
  const auto count = boost::sequential_vertex_coloring(lattice.graph(),
      boost::make_iterator_property_map(colors.begin(), get(boost::vertex_index, lattice.graph())));
  EXPECT_EQ(count, scenario.count);
  EXPECT_EQ(colors, scenario.colors);
  const auto edges = boost::edges(lattice.graph());
  for (auto edge = edges.first; edge != edges.second; ++edge) {
    const auto source = boost::source(*edge, lattice.graph());
    const auto target = boost::target(*edge, lattice.graph());
    EXPECT_NE(colors[source], colors[target]) << "bond " << source << " -> " << target;
  }
  for (const auto color : colors) EXPECT_LT(color, count);
}

// Preserve the four historical greedy-coloring results, including vertex order.
// The triangular result is four colors for this ordering, not a claim of optimality.
INSTANTIATE_TEST_SUITE_P(HistoricalLattices, LatticeColoring, ::testing::Values(
  ColoringScenario{"Chain", "chain lattice", 8, 2, {0, 1, 0, 1, 0, 1, 0, 1}},
  ColoringScenario{"Square", "square lattice", 4, 2, {0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0}},
  ColoringScenario{"FrustratedSquare", "frustrated square lattice", 4, 4, {0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3}},
  ColoringScenario{"Triangular", "triangular lattice", 4, 4, {0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3}}),
  [](const ::testing::TestParamInfo<ColoringScenario>& info) { return info.param.name; });
} // namespace
