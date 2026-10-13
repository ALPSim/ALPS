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

#include <alps/parapack/process.h>
#include "process_mpi_checks.hpp"
#include <stdexcept>
#include <vector>

TEST(ParallelProcessGroups, AllocatesExhaustsRecyclesAndReleasesGroups) {
  boost::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  const int group_size = world.size() >= 4 ? 2 : 1;
  const int group_count = world.size() / group_size;
  const int first_worker = world.size() % group_size;
  alps::process_helper_mpi process(world, group_size);
  if (world.rank() == 0) {
    EXPECT_EQ(process.num_groups(), group_count);
    EXPECT_EQ(process.num_free(), group_count);
    EXPECT_EQ(process.num_allocated(), 0);
    std::vector<alps::process_group> groups;
    for (int group = 0; group < group_count; ++group) {
      groups.push_back(process.allocate());
      EXPECT_EQ(groups.back().group_id, group);
      alps::ProcessList expected;
      for (int member = 0; member < group_size; ++member)
        expected.push_back(alps::Process(first_worker + group * group_size + member));
      EXPECT_EQ(groups.back().process_list, expected);
      EXPECT_EQ(int(groups.back().master()), first_worker + group * group_size);
      EXPECT_EQ(process.num_free(), group_count - group - 1);
      EXPECT_EQ(process.num_allocated(), group + 1);
    }
    EXPECT_THROW(process.allocate(), std::logic_error);
    process.release(groups.front());
    EXPECT_EQ(process.num_free(), 1);
    const auto recycled = process.allocate();
    EXPECT_EQ(recycled.group_id, groups.front().group_id);
    EXPECT_EQ(recycled.process_list, groups.front().process_list);
    EXPECT_EQ(process.num_free(), 0);
    for (std::size_t i = 1; i < groups.size(); ++i)
      process.release(groups[i]);
    process.release(recycled);
    EXPECT_EQ(process.num_free(), group_count);
    EXPECT_EQ(process.num_allocated(), 0);
    EXPECT_THROW(process.release(recycled), std::logic_error);
  }
  alps_test::expect_halted(process, world);
}
