/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010-2012 by Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/utility/vmusage.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <iostream>

#include <gtest/gtest.h>

TEST(VirtualMemory, QueryReportsNamedMetrics) {
  const auto metrics = alps::vmusage();
  // Platforms without /proc still expose the same keys with zero memory data.
  ASSERT_EQ(metrics.size(), 5u);
  for (const auto* name : {"Pid", "VmPeak", "VmSize", "VmHWM", "VmRSS"})
    ASSERT_EQ(metrics.count(name), 1u) << name;
  EXPECT_GE(metrics.at("VmPeak"), metrics.at("VmSize"));
  EXPECT_GE(metrics.at("VmHWM"), metrics.at("VmRSS"));
  const auto invalid = alps::vmusage(-2);
  ASSERT_EQ(invalid.size(), metrics.size());
  for (const auto& metric : invalid) EXPECT_EQ(metric.second, 0ul) << metric.first;
}
