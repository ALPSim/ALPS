/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/clone_timer.h>
#include <gtest/gtest.h>

TEST(CloneTimer, AdaptsBatchSizeWithoutWallClockSleeps) {
  using namespace boost::posix_time;
  alps::clone_timer timer(milliseconds(100));
  EXPECT_FALSE(timer.current_time().is_not_a_date_time());
  // Establish a synthetic timeline through the API's explicit time argument.
  const ptime initial(boost::gregorian::date(2000, 1, 1));
  timer.next_loops(1024, initial);
  EXPECT_EQ(timer.next_loops(1024, initial + milliseconds(10)), 2048);
  EXPECT_EQ(timer.next_loops(2048, initial + milliseconds(210)), 1024);
  EXPECT_EQ(timer.next_loops(1024, initial + milliseconds(310)), 1024);
  EXPECT_EQ(timer.next_loops(1, initial + seconds(1)), 1);
  timer.reset(.5);
  timer.next_loops(1024, initial);
  EXPECT_EQ(timer.next_loops(1024, initial + milliseconds(10)), 2048);
}
