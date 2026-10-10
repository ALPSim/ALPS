/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/queue.h>
#include <gtest/gtest.h>
#include <queue>

TEST(CheckQueue, OrdersScheduledChecksByEarliestTime) {
  using namespace boost::posix_time;
  const ptime initial(boost::gregorian::date(2000, 1, 1));
  EXPECT_LT(initial, initial + seconds(10));
  EXPECT_LT(initial + seconds(10), initial + minutes(5));
  EXPECT_LT(initial + minutes(5), initial + hours(1));
  std::priority_queue<alps::check_queue_element_t> queue;
  for (const int delay : {10, 2, 5, 1})
    queue.emplace(alps::check_type::checkpoint, initial + seconds(delay), 0, 0, 0);
  for (const int delay : {1, 2, 5, 10}) {
    ASSERT_FALSE(queue.empty());
    EXPECT_EQ(queue.top().time, initial + seconds(delay));
    queue.pop();
  }
  EXPECT_TRUE(queue.empty());
}
TEST(CheckQueue, CheckpointKeepsTaskCloneAndGroupIdentity) {
  using namespace boost::posix_time;
  const auto before = second_clock::local_time();
  const auto check = alps::next_checkpoint(2, 3, 4, seconds(10));
  const auto after = second_clock::local_time();
  EXPECT_EQ(check.type, alps::check_type::checkpoint);
  EXPECT_EQ(check.task_id, 2);
  EXPECT_EQ(check.clone_id, 3);
  EXPECT_EQ(check.group_id, 4);
  EXPECT_GE(check.time, before + seconds(10));
  EXPECT_LE(check.time, after + seconds(10));
}
