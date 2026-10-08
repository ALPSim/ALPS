/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/process.h>
#include "process_mpi_checks.hpp"

TEST(ParallelCommunicators, SplitsWorkersAndGroupHeads) {
  boost::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  ASSERT_LE(world.size(), 3);
  alps::process_helper_mpi process(world, 2);
  EXPECT_EQ(process.num_procs_per_group(), 2);
  EXPECT_EQ(process.num_total_processes(), world.size());
  EXPECT_EQ(process.comm_ctrl().size(), world.size());
  EXPECT_EQ(process.comm_ctrl().rank(), world.rank());

  // With three ranks the controller stays outside the two-worker group.
  const int first_worker = world.size() % 2;
  const bool worker = world.rank() >= first_worker;
  EXPECT_EQ(static_cast<bool>(process.comm_work()), worker);
  if (process.comm_work()) {
    EXPECT_EQ(process.comm_work().size(), 2);
    EXPECT_EQ(process.comm_work().rank(), world.rank() - first_worker);
  }
  const bool head = world.rank() == first_worker;
  EXPECT_EQ(static_cast<bool>(process.comm_head()), head);
  if (process.comm_head()) {
    EXPECT_EQ(process.comm_head().size(), 1);
    EXPECT_EQ(process.comm_head().rank(), 0);
  }
  if (world.rank() == 0) {
    EXPECT_EQ(process.num_groups(), 1);
    auto group = process.allocate();
    EXPECT_EQ(group.group_id, 0);
    EXPECT_EQ(group.process_list,
              (alps::ProcessList{alps::Process(first_worker), alps::Process(first_worker + 1)}));
    process.release(group);
  }
  alps_test::expect_halted(process, world);
}
