/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parameter/parameters.h>
#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

TEST(LegacyParametersHdf5, PreservesNumericPathsWhitespaceAndUnresolvedValues) {
  std::istringstream input(R"PARAMS(a = 2
b = "/home/ALPS/lib/xml"
c = 'abc def '
d = "../lattices.xml"
f= .1.4
g= ....5
)PARAMS");
  const alps::Parameters original(input);
  const std::vector<std::pair<std::string, std::string>> expected{
      {"a", "2"},        {"b", "/home/ALPS/lib/xml"},
      {"c", "abc def "}, {"d", "../lattices.xml"},
      {"f", ".1.4"},     {"g", "....5"}};
  alps::testing::TemporaryDirectory directory;
  const auto filename = (directory.path() / "parameters.h5").string();
  {
    alps::hdf5::archive archive(filename, "w");
    archive["/parameters"] << original;
  }
  alps::Parameters restored;
  {
    alps::hdf5::archive archive(filename, "r");
    archive["/parameters"] >> restored;
  }
  ASSERT_EQ(restored.size(), expected.size());
  for (const auto &entry : expected) {
    SCOPED_TRACE(entry.first);
    ASSERT_TRUE(restored.defined(entry.first));
    EXPECT_EQ(std::string(original[entry.first].c_str()), entry.second);
    EXPECT_EQ(std::string(restored[entry.first].c_str()), entry.second);
  }
}
