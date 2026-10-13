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

#include <alps/parameter/parameter.h>
#include <gtest/gtest.h>
#include <string>

struct Assignment {
  const char *input;
  const char *key;
  const char *value;
};
void PrintTo(const Assignment &sample, std::ostream *output) { *output << sample.key; }
class LegacyParameter : public ::testing::TestWithParam<Assignment> {};
TEST_P(LegacyParameter, ParsesKeysValuesQuotesAndEnvironment) {
  const auto sample = GetParam();
  SCOPED_TRACE(sample.input);
  alps::Parameter parameter;
  ASSERT_NO_THROW(parameter.parse(sample.input));
  EXPECT_EQ(parameter.key(), sample.key);
  EXPECT_EQ(std::string(parameter.value().c_str()), sample.value);
}
INSTANTIATE_TEST_SUITE_P(
    HistoricalInputs, LegacyParameter,
    ::testing::Values(Assignment{"L=10 ; ", "L", "10"}, Assignment{"M=1 ", "M", "1"},
                      Assignment{"N =", "N", ""}, Assignment{"O = ;", "O", ""},
                      Assignment{"T=0.1 ;", "T", "0.1"}, Assignment{"beta='1/T';", "beta", "1/T"},
                      Assignment{"a0 = \"abc def \";", "a0", "abc def "},
                      Assignment{"a# = 'abc def '", "a#", "abc def "},
                      Assignment{"a[test 2] = 3;", "a[test 2]", "3"},
                      Assignment{"a[ test 3 ] = 3;", "a[ test 3 ]", "3"},
                      Assignment{"dir = \"${DIR}/lib/xml${ DIR}\"", "dir",
                                 "/home/alps/lib/xml${ DIR}"},
                      Assignment{"b = sqrt(4)+2*Pi/L-X^5", "b", "sqrt(4)+2*Pi/L-X^5"},
                      Assignment{"c = 'a - b[4]'", "c", "a - b[4]"}));
