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

#include <alps/parapack/clone_info.h>
#include <boost/mpi.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <vector>

TEST(ParallelCloneInfo, KeepsRankSpecificSeedsCheckpointsAndMasterPhases) {
  boost::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  alps::Parameters params;
  params["SEED"] = 29832;
  alps::clone_info_mpi info(world, 0, params, "info_test");
  alps::clone_info_mpi same(world, 0, params, "same");
  alps::clone_info_mpi other(world, 1, params, "other");
  EXPECT_EQ(info.clone_id(), 0);
  EXPECT_TRUE(info.has_seed());
  EXPECT_EQ(info.worker_seed(), same.worker_seed());
  EXPECT_NE(info.worker_seed(), other.worker_seed());
  std::vector<alps::seed_t> workers, disorder;
  boost::mpi::all_gather(world, info.worker_seed(), workers);
  boost::mpi::all_gather(world, info.disorder_seed(), disorder);
  std::sort(workers.begin(), workers.end());
  EXPECT_EQ(std::adjacent_find(workers.begin(), workers.end()), workers.end());
  EXPECT_TRUE(std::all_of(disorder.begin(), disorder.end(),
                          [&](alps::seed_t seed) { return seed == info.disorder_seed(); }));
  EXPECT_EQ(info.dumpfile(), "info_test.clone1.worker" + std::to_string(world.rank() + 1));
  EXPECT_EQ(info.checkpoints().size(), world.rank() == 0 ? world.size() : 1u);
  if (world.rank() == 0) {
    for (int rank = 0; rank < world.size(); ++rank) {
      EXPECT_EQ(info.checkpoints()[rank], "info_test.clone1.worker" + std::to_string(rank + 1));
    }
  }
  info.start("test 1");
  info.stop();
  info.start("test 2");
  info.stop();
  if (world.rank() == 0) {
    ASSERT_EQ(info.phases().size(), 2);
    EXPECT_EQ(info.phases()[0].phase(), "test 1");
    EXPECT_EQ(info.phases()[1].phase(), "test 2");
    EXPECT_EQ(info.hosts().size(), world.size());
    for (const auto &host : info.hosts())
      EXPECT_FALSE(host.empty());
    EXPECT_GE(info.elapsed().total_microseconds(), 0);
  } else {
    EXPECT_TRUE(info.phases().empty());
  }
  EXPECT_DOUBLE_EQ(info.progress(), 0);
  info.set_progress(.593483);
  EXPECT_DOUBLE_EQ(info.progress(), .593483);
  info.set_progress(1);
  EXPECT_DOUBLE_EQ(info.progress(), 1);
}
