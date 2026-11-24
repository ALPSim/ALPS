// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The cluster library reader (alps_graph::lookup_basis_from_name, #443). It
// matched dca.LATTICE as a substring of the first "... superlattice" line, so
// "4A" selected Betts2D-4A, and it read the basis entries into ints without
// checking the stream, so "2.5" became 2 and "x" or a missing entry was left
// uninitialized.

#include "gtest/gtest.h"
#include "cluster_graph.h"

#include <cstdio>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<std::string> superlattice_names(const std::string& library) {
  std::ifstream in(library);
  std::stringstream text;
  text << in.rdbuf();
  const std::string s = text.str();
  const std::regex name("<LATTICE\\s+name\\s*=\\s*\"([^\"]*) superlattice\"");
  std::vector<std::string> names;
  for (auto it = std::sregex_iterator(s.begin(), s.end(), name); it != std::sregex_iterator(); ++it)
    names.push_back((*it)[1]);
  return names;
}

// The number of sites in names like Betts2D-8A, Kent3D-16Z or Chain1D_7, 0 otherwise
int sites_from_name(const std::string& name) {
  std::smatch m;
  if (std::regex_match(name, m, std::regex("[A-Za-z]+[123]D[-_]([0-9]+)[A-Z]?")))
    return std::stoi(m[1]);
  return 0;
}

std::string error_of(const std::string& name, const std::string& library) {
  try {
    alps_graph g(name, library, 2);
  } catch (std::exception const& e) {
    return e.what();
  }
  return "";
}

void expect_contains(const std::string& text, const std::string& part) {
  EXPECT_NE(text.find(part), std::string::npos) << "expected \"" << part << "\" in: " << text;
}

const std::string dir = DCA_CLUSTER_DIR;
}

// Every cluster of the shipped libraries is read, with the dimension and the
// number of sites its name says. "chain" in cluster.xml is a parametrized
// lattice of the ALPS lattice library and has no numeric basis.
TEST(ClusterLibrary, ReadsEveryShippedCluster) {
  for (const char* file : {"cluster.xml", "cluster1d.xml", "cluster3d.xml"}) {
    const std::string library = dir + "/" + file;
    const std::vector<std::string> names = superlattice_names(library);
    EXPECT_GT(names.size(), 10u) << library;
    int with_size = 0;
    for (const std::string& name : names) {
      if (name == "chain") continue;
      SCOPED_TRACE(name + " in " + file);
      alps_graph g(name, library, 2);
      const int sites = sites_from_name(name);
      if (sites) {
        ++with_size;
        EXPECT_EQ(g.unit_cell_volume(), sites);
        EXPECT_EQ(g.dimension(), name[name.find('D') - 1] - '0');
      }
    }
    EXPECT_GT(with_size, 10) << library;
  }
  alps_graph dimer("dimer", dir + "/cluster.xml", 2);
  EXPECT_EQ(dimer.unit_cell_volume(), 2);
  EXPECT_EQ(dimer.dimension(), 2);
}

TEST(ClusterLibrary, NameMustMatchExactly) {
  const std::string library = dir + "/cluster.xml";
  for (const char* partial : {"4A", "D-8A", "Betts2D-8", "betts2d-8a", "Betts2D-8A "})
    expect_contains(error_of(partial, library),
                    std::string("cluster ") + partial + " in " + library + ": not found. dca.LATTICE has to be the name");
  expect_contains(error_of("chain", library), "basis vector \"SITES\": \"SITES\" is not an integer.");
  expect_contains(error_of("chain", library), "Chain1D_<N> in cluster1d.xml");
}

TEST(ClusterLibrary, RefusesMalformedEntries) {
  const std::string library = "cluster_library_test.xml";
  std::ofstream(library) <<
      "<LATTICES>\n"
      "<LATTICE name=\"letter superlattice\" dimension=\"2\"><BASIS><VECTOR>2 x</VECTOR><VECTOR>0 2</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"half superlattice\" dimension=\"2\"><BASIS><VECTOR>2.5 0</VECTOR><VECTOR>0 2</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"short superlattice\" dimension=\"2\"><BASIS><VECTOR>2</VECTOR><VECTOR>0 2</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"one superlattice\" dimension=\"2\"><BASIS><VECTOR>2 0</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"flat superlattice\" dimension=\"2\"><BASIS><VECTOR>2 2</VECTOR><VECTOR>1 1</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"four superlattice\" dimension=\"4\"><BASIS><VECTOR>1 0 0 0</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"twice superlattice\" dimension=\"1\"><BASIS><VECTOR>2</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"twice superlattice\" dimension=\"1\"><BASIS><VECTOR>3</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"nobasis superlattice\" dimension=\"1\"></LATTICE>\n"
      // the layouts the old parser required are not needed: attributes in
      // another order, spaces around '=', the basis on one line, 2.0
      "<LATTICE dimension = \"2\"\n   name = \"layout superlattice\">\n<BASIS><VECTOR> 2.0  2 </VECTOR>\n"
      "<VECTOR>2 -2</VECTOR></BASIS></LATTICE>\n"
      "<LATTICE name=\"left superlattice\" dimension=\"1\"><BASIS><VECTOR>-3</VECTOR></BASIS></LATTICE>\n"
      // self-closing tag: the next lattice's basis must not be used
      "<LATTICE name=\"selfclose superlattice\" dimension=\"2\"/>\n"
      "<LATTICE name=\"afterself superlattice\" dimension=\"2\"><BASIS><VECTOR>3 0</VECTOR><VECTOR>0 3</VECTOR></BASIS></LATTICE>\n"
      // commented-out entry: must not be found, and a commented duplicate must
      // not count towards "defines it more than once"
      "<LATTICE name=\"single superlattice\" dimension=\"1\"><BASIS><VECTOR>5</VECTOR></BASIS></LATTICE>\n"
      "<!-- <LATTICE name=\"incomment superlattice\" dimension=\"1\"><BASIS><VECTOR>2</VECTOR></BASIS></LATTICE> -->\n"
      "<!-- <LATTICE name=\"single superlattice\" dimension=\"1\"><BASIS><VECTOR>7</VECTOR></BASIS></LATTICE> -->\n"
      // overflow: 1e300 rounds to itself but overflows int
      "<LATTICE name=\"huge superlattice\" dimension=\"1\"><BASIS><VECTOR>1e300</VECTOR></BASIS></LATTICE>\n"
      "</LATTICES>\n";
  const std::string where = " in " + library + ": ";
  expect_contains(error_of("letter", library), "cluster letter" + where + "basis vector \"2 x\": \"x\" is not an integer.");
  expect_contains(error_of("half", library), "basis vector \"2.5 0\": \"2.5\" is not an integer.");
  expect_contains(error_of("short", library), "basis vector \"2\" has 1 entries, dimension=\"2\" needs 2.");
  expect_contains(error_of("one", library), "the basis has 1 vectors, dimension=\"2\" needs 2.");
  expect_contains(error_of("flat", library), "the basis vectors are linearly dependent");
  expect_contains(error_of("four", library), "dimension=\"4\"; clusters have dimension 1, 2 or 3.");
  expect_contains(error_of("twice", library), "the library defines it more than once.");
  expect_contains(error_of("nobasis", library), "has no <BASIS> ... </BASIS>.");
  expect_contains(error_of("selfclose", library), "has no <BASIS> ... </BASIS>.");
  expect_contains(error_of("incomment", library), "not found. dca.LATTICE has to be the name");
  // a commented-out duplicate must not trigger "defines it more than once"
  EXPECT_EQ(alps_graph("single", library, 1).unit_cell_volume(), 5);
  expect_contains(error_of("huge", library), "is not an integer.");

  alps_graph layout("layout", library, 2);
  EXPECT_EQ(layout.unit_cell_volume(), 8);
  ASSERT_EQ(layout.basis_vectors().size(), 2u);
  EXPECT_EQ(layout.basis_vectors()[0], (std::vector<double>{2, 2}));
  EXPECT_EQ(layout.basis_vectors()[1], (std::vector<double>{2, -2}));
  // the volume of a 1D cluster is the length of its basis vector
  EXPECT_EQ(alps_graph("left", library, 2).unit_cell_volume(), 3);
  std::remove(library.c_str());
}
