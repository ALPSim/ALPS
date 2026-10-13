#include <gtest/gtest.h>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/config.h>
#include <cstddef>
#include <iostream>

TEST(PlatformTypes, FixedWidthIntegersMatchWireSizes) {
  EXPECT_EQ(sizeof(alps::int8_t), 1u);
  EXPECT_EQ(sizeof(alps::int16_t), 2u);
  EXPECT_EQ(sizeof(alps::int32_t), 4u);
  EXPECT_EQ(sizeof(alps::int64_t), 8u);
  EXPECT_GE(sizeof(std::size_t), 4u);
}
