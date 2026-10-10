#include <gtest/gtest.h>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/utility/os.hpp>
#include <alps/version.h>
#include <iostream>

TEST(OperatingSystem, QueriesAreAvailable) {
  EXPECT_FALSE(alps::hostname().empty());
  EXPECT_FALSE(alps::username().empty());
}
