/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/parapack.h>
#include <gtest/gtest.h>
#include <regex>

TEST(ParapackVersion, ReportsANumericRelease) {
  EXPECT_TRUE(std::regex_search(alps::parapack::alps_version(), std::regex("[0-9]+\\.[0-9]+")));
}
