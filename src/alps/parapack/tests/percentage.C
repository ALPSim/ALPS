/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/util.h>
#include <gtest/gtest.h>

TEST(Percentage, AcceptsWhitespaceAroundPercentSign) {
  for (const auto *input : {"10%", "10 %", " 10%", " 10 %  "}) {
    SCOPED_TRACE(input);
    EXPECT_DOUBLE_EQ(alps::parse_percentage(input), .1);
  }
  EXPECT_DOUBLE_EQ(alps::parse_percentage("0.1%"), .001);
}
