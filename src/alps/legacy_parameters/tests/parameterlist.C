/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parameter/parameterlist.h>
#include <alps/osiris/xdrdump.h>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
using ParameterMap = std::map<std::string, std::string>;
const std::vector<ParameterMap> expected{
    {{"A", "0.1"}, {"B", "0.2"}, {"C", "200"}},
    {{"A", "0.1"}, {"D", "0.2"}, {"dir", "/home/alps/lib/xml"}, {"E", "0.2"}, {"F", "100"}},
    {{"A", "0.1"}, {"D", "0.2"}, {"dir", "/home/alps/lib/xml"}, {"G", "0.5"}, {"H", "100"}},
    {{"A", "0.1"}, {"D", "0.2"}, {"dir", "/home/alps/lib/xml"}},
    {{"G", "0.4"}, {"H", "200"}},
    {{"H", "300"}}};
void expect_parameters(const alps::ParameterList &params) {
  ASSERT_EQ(params.size(), expected.size());
  for (std::size_t i = 0; i < params.size(); ++i) {
    SCOPED_TRACE(i);
    ParameterMap values;
    for (const auto &value : params[i])
      values.emplace(value.key(), value.value().c_str());
    EXPECT_EQ(values, expected[i]);
  }
}
} // namespace
TEST(LegacyParameterList, GlobalInheritanceClearStopAndXdrRoundTrip) {
  std::istringstream input(R"PARAMS(/* C-style comment */

// C++-style comments are also accepted

A=0.1;
{ B=0.2; C=200; }

D=0.2;
dir = "${DIR}/lib/xml"

{

E=0.2, F=100
}

{ G=0.5; H=100}

{}

#clear // clear global parameters

{ G=0.4; H=200 }

{ /* G=0.2; */ H=300; }

#stop // stop reading here

{ G=0.3; H=400; }
)PARAMS");
  alps::ParameterList params(input);
  expect_parameters(params);
  alps::testing::TemporaryDirectory directory;
  const boost::filesystem::path dump((directory.path() / "parameters.xdr").string());
  {
    alps::OXDRFileDump output(dump);
    output << params;
  }
  params.clear();
  {
    alps::IXDRFileDump input(dump);
    input >> params;
  }
  expect_parameters(params);
}
