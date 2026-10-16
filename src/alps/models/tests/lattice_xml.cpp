// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The XML descriptions of kit lattices: the <GRAPH> and <LATTICES> documents,
// content-addressed custom names and the catalog graph name a boundary
// selects. The string checks follow skilledwolf/alpsx (branch
// agent/v3-execution-and-api-audit at c99b526ec^, tests/alps/test_lattice.py
// TestGraphXML); the round trip through the legacy reader is in the legacy
// model module's lattice_xml_bridge test.
#include <alps/models/lattice.hpp>
#include <alps/models/lattice_xml.hpp>
#include <gtest/gtest.h>

#include <regex>
#include <string>

namespace lat = alps::lattice;
using lat::Boundary;

namespace {
bool contains(std::string const& text, std::string const& fragment) { return text.find(fragment) != std::string::npos; }
} // namespace

TEST(LatticeXml, GraphDocumentListsSitesWithCoordinatesAndBonds) {
    auto const chain = lat::to_graph_xml(lat::chain(4, 1.0, Boundary::open));
    EXPECT_TRUE(contains(chain, "<GRAPH vertices=\"4\" edges=\"3\">"));
    EXPECT_TRUE(contains(chain, "<VERTEX id=\"1\" type=\"0\"><COORDINATE>0</COORDINATE></VERTEX>"));
    EXPECT_TRUE(contains(chain, "<VERTEX id=\"4\" type=\"0\"><COORDINATE>3</COORDINATE></VERTEX>"));
    EXPECT_TRUE(contains(chain, "<EDGE source=\"1\" target=\"2\" type=\"0\"/>"));
    EXPECT_TRUE(contains(chain, "<EDGE source=\"3\" target=\"4\" type=\"0\"/>"));
    EXPECT_FALSE(contains(chain, "<EDGE source=\"4\" target=\"1\""));
    EXPECT_TRUE(contains(lat::to_graph_xml(lat::square(2, 2, 1.0, Boundary::open)), "<GRAPH vertices=\"4\" edges=\"4\">"));
    auto const honeycomb = lat::to_graph_xml(lat::honeycomb(2, 2, 1.0, Boundary::open));
    EXPECT_TRUE(contains(honeycomb, "type=\"0\"><COORDINATE>"));
    EXPECT_TRUE(contains(honeycomb, "type=\"1\"><COORDINATE>"));
    auto const ladder = lat::to_graph_xml(lat::ladder(3, 1.0, Boundary::open));
    EXPECT_TRUE(contains(ladder, "<EDGE source=\"1\" target=\"3\" type=\"0\"/>")); // leg
    EXPECT_TRUE(contains(ladder, "<EDGE source=\"1\" target=\"2\" type=\"1\"/>")); // rung
    // Real-space coordinates, 17 significant digits.
    auto const triangular = lat::to_graph_xml(lat::triangular(2, 2));
    EXPECT_TRUE(contains(triangular, "<COORDINATE>0.5 0.8660254037844386</COORDINATE>"));
}

TEST(LatticeXml, LatticeGraphDocumentBakesCellCoordinatesExtentAndBoundary) {
    auto const xml = lat::to_latticegraph_xml(lat::honeycomb(2, 3), "x");
    EXPECT_TRUE(contains(xml, "<LATTICES>\n<LATTICE name=\"x base\" dimension=\"2\"><BASIS>"
                              "<VECTOR>1 0</VECTOR><VECTOR>0.5 0.8660254037844386</VECTOR></BASIS></LATTICE>\n"));
    EXPECT_TRUE(contains(xml, "<UNITCELL name=\"x cell\" dimension=\"2\">\n"
                              "  <VERTEX type=\"0\"><COORDINATE>0 0</COORDINATE></VERTEX>\n"
                              "  <VERTEX type=\"1\"><COORDINATE>-0.33333333333333331 0.66666666666666663</COORDINATE></VERTEX>\n"
                              "  <EDGE type=\"0\"><SOURCE vertex=\"1\" offset=\"0 0\"/><TARGET vertex=\"2\" offset=\"0 0\"/></EDGE>\n"
                              "  <EDGE type=\"1\"><SOURCE vertex=\"2\" offset=\"0 0\"/><TARGET vertex=\"1\" offset=\"0 1\"/></EDGE>\n"
                              "  <EDGE type=\"2\"><SOURCE vertex=\"1\" offset=\"0 0\"/><TARGET vertex=\"2\" offset=\"1 -1\"/></EDGE>\n"
                              "</UNITCELL>\n"));
    EXPECT_TRUE(contains(xml, "<LATTICEGRAPH name=\"x\">\n  <FINITELATTICE>\n    <LATTICE ref=\"x base\"/>\n"
                              "    <EXTENT dimension=\"1\" size=\"2\"/>\n    <EXTENT dimension=\"2\" size=\"3\"/>\n"
                              "    <BOUNDARY dimension=\"1\" type=\"periodic\"/>\n    <BOUNDARY dimension=\"2\" type=\"periodic\"/>\n"
                              "  </FINITELATTICE>\n  <UNITCELL ref=\"x cell\"/>\n</LATTICEGRAPH>\n</LATTICES>\n"));
    lat::UnitCell cell(2);
    cell.add_vertex({0., 0.});
    cell.add_edge(0, 0, 0, {}, {1, 0});
    cell.add_edge(0, 0, 0, {}, {0, 1});
    auto const cylinder = lat::to_latticegraph_xml(cell.tile({4, 3}, {}, {Boundary::periodic, Boundary::open}));
    EXPECT_TRUE(contains(cylinder, "<LATTICEGRAPH name=\"alps_custom_lattice\">"));
    EXPECT_TRUE(contains(cylinder, "<BOUNDARY dimension=\"1\" type=\"periodic\"/>\n    <BOUNDARY dimension=\"2\" type=\"open\"/>"));
    auto const cubic = lat::to_latticegraph_xml(lat::simple_cubic(2, 3, 4, 1.0, Boundary::open), "c");
    EXPECT_TRUE(contains(cubic, "<EXTENT dimension=\"3\" size=\"4\"/>"));
    EXPECT_TRUE(contains(cubic, "<BOUNDARY dimension=\"3\" type=\"open\"/>"));
}

TEST(LatticeXml, CustomNamesAreContentAddressed) {
    std::regex const shape("alps_custom_lattice_[0-9a-f]{16}");
    auto const a = lat::custom_lattice_name(lat::honeycomb(2, 3));
    EXPECT_TRUE(std::regex_match(a, shape)) << a;
    EXPECT_EQ(a, lat::custom_lattice_name(lat::honeycomb(2, 3)));
    EXPECT_NE(a, lat::custom_lattice_name(lat::honeycomb(3, 2)));
    EXPECT_NE(a, lat::custom_lattice_name(lat::honeycomb(2, 3, 1.0, Boundary::open)));
    EXPECT_NE(a, lat::custom_lattice_name(lat::honeycomb(2, 3, 2.0)));
    EXPECT_NE(a, lat::custom_lattice_name(lat::triangular(2, 3)));
    // The name does not depend on the name the document is emitted under.
    EXPECT_NE(lat::to_latticegraph_xml(lat::honeycomb(2, 3), "p"), lat::to_latticegraph_xml(lat::honeycomb(2, 3), "q"));
    EXPECT_TRUE(contains(lat::to_latticegraph_xml(lat::honeycomb(2, 3), a), "<LATTICEGRAPH name=\"" + a + "\">"));
}

TEST(LatticeXml, CatalogGraphNamesFollowTheBoundary) {
    EXPECT_EQ(lat::catalog_graph_name(lat::chain(4)), "chain lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::chain(4, 1.0, Boundary::open)), "open chain lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::nnn_chain(4, 1.0, Boundary::open)), "nnn open chain lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::square(3, 3)), "square lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::square(3, 3, 1.0, Boundary::open)), "open square lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::simple_cubic(2, 2, 2, 1.0, Boundary::open)), "open simple cubic lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::ladder(4)), "ladder");
    EXPECT_EQ(lat::catalog_graph_name(lat::ladder(4, 1.0, Boundary::open)), "open ladder");
    EXPECT_EQ(lat::catalog_graph_name(lat::honeycomb(2, 2)), "honeycomb lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::kagome(2, 2)), "kagome lattice");
    EXPECT_EQ(lat::catalog_graph_name(lat::triangular(3, 3, 1.0, Boundary::open)), "");
    EXPECT_EQ(lat::catalog_graph_name(lat::honeycomb(2, 2, 1.0, Boundary::open)), "");
    lat::UnitCell cell(2);
    cell.add_vertex({0., 0.});
    cell.add_edge(0, 0, 0, {}, {1, 0});
    cell.add_edge(0, 0, 0, {}, {0, 1});
    EXPECT_EQ(lat::catalog_graph_name(cell.tile({4, 3})), "");
    auto const cylinder = lat::Lattice::build(cell, {{1., 0.}, {0., 1.}}, {4, 3},
                                              {Boundary::periodic, Boundary::open}, "square lattice");
    EXPECT_EQ(lat::catalog_graph_name(cylinder), "");
}
