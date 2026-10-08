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
  alps::process_helper_mpi process(world, 1);
  if (world.rank() == 0) {
    EXPECT_EQ(process.num_groups(), world.size());
    EXPECT_EQ(process.num_free(), world.size());
    EXPECT_EQ(process.num_allocated(), 0);
    std::vector<alps::process_group> groups;
    for (int rank = 0; rank < world.size(); ++rank) {
      groups.push_back(process.allocate());
      EXPECT_EQ(groups.back().group_id, rank);
      EXPECT_EQ(groups.back().process_list, alps::ProcessList{alps::Process(rank)});
      EXPECT_EQ(int(groups.back().master()), rank);
      EXPECT_EQ(process.num_free(), world.size() - rank - 1);
      EXPECT_EQ(process.num_allocated(), rank + 1);
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
    EXPECT_EQ(process.num_free(), world.size());
    EXPECT_EQ(process.num_allocated(), 0);
    EXPECT_THROW(process.release(recycled), std::logic_error);
  }
  alps_test::expect_halted(process, world);
}
