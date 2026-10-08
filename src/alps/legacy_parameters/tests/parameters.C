/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parameter/parameters.h>
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <vector>

TEST(LegacyParameters, ParsesCommentsPreservesOrderAndSupportsMutation) {
  std::istringstream input(R"PARAMS(L=10; M= 1,

T=0.1;
beta="1/T"; // comment
dir = "${DIR}/lib/xml", a0 = "abc def ";
a0 = 'abc def '
a[test 2] = 3, /* comment */  a[ test 3 ] = 3;
b = sqrt(4)+2*Pi/L-X^5;


)PARAMS");
  alps::Parameters params(input);
  const std::vector<std::pair<std::string, std::string>> expected{{"L", "10"},
                                                                  {"M", "1"},
                                                                  {"T", "0.1"},
                                                                  {"beta", "1/T"},
                                                                  {"dir", "/home/alps/lib/xml"},
                                                                  {"a0", "abc def "},
                                                                  {"a[test 2]", "3"},
                                                                  {"a[ test 3 ]", "3"},
                                                                  {"b", "sqrt(4)+2*Pi/L-X^5"}};
  ASSERT_EQ(params.size(), expected.size());
  auto actual = params.begin();
  for (const auto &entry : expected) {
    EXPECT_EQ(actual->key(), entry.first);
    EXPECT_EQ(std::string(actual->value().c_str()), entry.second);
    ++actual;
  }
  params["L"] = 3;
  EXPECT_EQ(std::string(params["L"].c_str()), "3");
  params.erase("a0");
  EXPECT_FALSE(params.defined("a0"));
  EXPECT_EQ(params.size(), expected.size() - 1);
  alps::Parameters copy(params);
  copy["N"] = copy["L"];
  EXPECT_EQ(std::string(copy["N"].c_str()), "3");
  EXPECT_FALSE(params.defined("N"));
}
