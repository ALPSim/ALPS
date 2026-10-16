/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/lattice.h>
#include <sstream>
#include <string>

namespace {
std::string round_trip(const std::string& input) {
  std::istringstream in(input);
  alps::XMLTag tag = alps::parse_tag(in);
  alps::InhomogeneityDescriptor descriptor(tag, in);
  std::ostringstream out;
  {
    alps::oxstream xml(out);
    descriptor.write_xml(xml);
  }
  return out.str();
}
} // namespace

// Regression test for ALPSim/ALPS#187: the TARGET of a changed edge was
// stored in source_, losing the SOURCE and leaving target_ uninitialised.
TEST(LatticeDisorder, ChangedEdgeKeepsSourceAndTarget) {
  const std::string output = round_trip(
      "<CHANGED><EDGE type=\"2\">"
      "<SOURCE cell=\"1 0\" vertex=\"3\"/>"
      "<TARGET cell=\"0 1\" offset=\"1 1\" vertex=\"7\"/>"
      "</EDGE></CHANGED></LATTICEGRAPH>");

  const std::string::size_type source = output.find("<SOURCE");
  const std::string::size_type target = output.find("<TARGET");
  ASSERT_NE(source, std::string::npos) << output;
  ASSERT_NE(target, std::string::npos) << output;

  const std::string source_xml = output.substr(source, target - source);
  const std::string target_xml = output.substr(target);
  EXPECT_NE(source_xml.find("cell=\"1 0\""), std::string::npos) << output;
  EXPECT_NE(source_xml.find("vertex=\"3\""), std::string::npos) << output;
  EXPECT_EQ(source_xml.find("offset="), std::string::npos) << output;
  EXPECT_NE(target_xml.find("cell=\"0 1\""), std::string::npos) << output;
  EXPECT_NE(target_xml.find("offset=\"1 1\""), std::string::npos) << output;
  EXPECT_NE(target_xml.find("vertex=\"7\""), std::string::npos) << output;
}

// A changed vertex must keep its cell, offset and vertex attributes, and be
// written back as <VERTEX> so the output can be read again.
TEST(LatticeDisorder, ChangedVertexRoundTrips) {
  const std::string input =
      "<CHANGED><VERTEX type=\"4\" cell=\"2 1\" offset=\"0 1\" vertex=\"5\"/>"
      "</CHANGED></LATTICEGRAPH>";
  const std::string output = round_trip(input);

  EXPECT_NE(output.find("<VERTEX"), std::string::npos) << output;
  EXPECT_EQ(output.find("<CELL"), std::string::npos) << output;
  EXPECT_NE(output.find("type=\"4\""), std::string::npos) << output;
  EXPECT_NE(output.find("cell=\"2 1\""), std::string::npos) << output;
  EXPECT_NE(output.find("offset=\"0 1\""), std::string::npos) << output;
  EXPECT_NE(output.find("vertex=\"5\""), std::string::npos) << output;

  EXPECT_EQ(round_trip(output + "</LATTICEGRAPH>"), output);
}

TEST(LatticeDisorder, ChangedVertexRequiresCell) {
  EXPECT_THROW(round_trip("<CHANGED><VERTEX type=\"4\" vertex=\"5\"/>"
                          "</CHANGED></LATTICEGRAPH>"),
               std::runtime_error);
}

TEST(LatticeDisorder, DefaultVertexReferenceIsZeroed) {
  const alps::detail::BasicVertexReference reference;
  EXPECT_EQ(reference.vertex(), 0u);
  EXPECT_TRUE(reference.cell_offset().empty());
  EXPECT_TRUE(reference.offset().empty());
}
