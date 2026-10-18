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

TEST(LegacyParametersHdf5, StoresUnevaluableValuesAsWritten) {
  std::istringstream input(R"PARAMS(SOLVER = "Hirsch-Fye"
U_MATRIX = "u-matrix.dat"
SHIFT = "U/2+JH"
PHASE = "x*cos(2*pi)"
EMPTY = ""
BLANK = "  "
U = 3
)PARAMS");
  const alps::Parameters original(input);
  alps::testing::TemporaryDirectory directory;
  const auto filename = (directory.path() / "parameters.h5").string();
  {
    alps::hdf5::archive archive(filename, "w");
    archive["/parameters"] << original;
  }
  alps::hdf5::archive archive(filename, "r");
  for (const char *key : {"SOLVER", "U_MATRIX", "SHIFT", "PHASE", "EMPTY", "BLANK"}) {
    SCOPED_TRACE(key);
    const std::string path = std::string("/parameters/") + key;
    ASSERT_TRUE(archive.is_datatype<std::string>(path));
    std::string stored;
    archive[path] >> stored;
    EXPECT_EQ(stored, std::string(original[key].c_str()));
  }
}

TEST(LegacyParametersHdf5, StoresEvaluableValuesAsNumbers) {
  std::istringstream input(R"PARAMS(L = 16
NEGATIVE = -3
N = 2*L
t = 0.707106781186547
THIRD = 1/3
MU = U/2
U = 3
INT_MAX = 2147483647
SEED = 4294967295
HUGE = -1e300
)PARAMS");
  const alps::Parameters original(input);
  alps::testing::TemporaryDirectory directory;
  const auto filename = (directory.path() / "parameters.h5").string();
  {
    alps::hdf5::archive archive(filename, "w");
    archive["/parameters"] << original;
  }
  alps::hdf5::archive archive(filename, "r");
  const std::vector<std::pair<std::string, int>> integers{
      {"L", 16}, {"NEGATIVE", -3}, {"N", 32}, {"U", 3}, {"INT_MAX", 2147483647}};
  for (const auto &entry : integers) {
    SCOPED_TRACE(entry.first);
    ASSERT_TRUE(archive.is_datatype<int>("/parameters/" + entry.first));
    int stored = 0;
    archive["/parameters/" + entry.first] >> stored;
    EXPECT_EQ(stored, entry.second);
  }
  // Integral values outside the int range keep full precision as doubles.
  const std::vector<std::pair<std::string, double>> reals{
      {"t", 0.707106781186547}, {"THIRD", 1. / 3.}, {"MU", 1.5},
      {"SEED", 4294967295.}, {"HUGE", -1e300}};
  for (const auto &entry : reals) {
    SCOPED_TRACE(entry.first);
    ASSERT_TRUE(archive.is_datatype<double>("/parameters/" + entry.first));
    double stored = 0;
    archive["/parameters/" + entry.first] >> stored;
    EXPECT_EQ(stored, entry.second);
  }
}
