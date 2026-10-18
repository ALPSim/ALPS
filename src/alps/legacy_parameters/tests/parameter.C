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
#include <sstream>
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

struct Written {
  std::string value;
  std::string text;
};
void PrintTo(const Written &sample, std::ostream *output) { *output << sample.text; }
class LegacyParameterOutput : public ::testing::TestWithParam<Written> {};

// Values are quoted exactly when the unquoted form would not parse back.
TEST_P(LegacyParameterOutput, WritesTextThatParsesBack) {
  const auto sample = GetParam();
  std::ostringstream output;
  output << alps::Parameter("KEY", sample.value);
  EXPECT_EQ(output.str(), "KEY = " + sample.text + ";");
  alps::Parameter parsed;
  ASSERT_NO_THROW(parsed.parse(output.str(), false));
  EXPECT_EQ(parsed.key(), "KEY");
  EXPECT_EQ(std::string(parsed.value().c_str()), sample.value);
}
INSTANTIATE_TEST_SUITE_P(
    Values, LegacyParameterOutput,
    ::testing::Values(Written{"16", "16"}, Written{"-2.5E+2", "-2.5E+2"},
                      Written{"Hirsch-Fye", "Hirsch-Fye"}, Written{"U/2+JH", "U/2+JH"},
                      Written{"sqrt(4)*Pi^2", "sqrt(4)*Pi^2"},
                      Written{"../g0/tau-input.dat", "../g0/tau-input.dat"},
                      Written{"it's", "it's"}, Written{"chain lattice", "\"chain lattice\""},
                      Written{" padded ", "\" padded \""}, Written{"Sz:Sz", "\"Sz:Sz\""},
                      Written{"up,down", "\"up,down\""}, Written{"a;b", "\"a;b\""},
                      Written{"x=y", "\"x=y\""}, Written{"a[1]", "\"a[1]\""},
                      Written{"${HOME}", "\"${HOME}\""}, Written{"", "\"\""},
                      Written{"say \"hi\"", "'say \"hi\"'"}));
