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

TEST(ParallelHalt, WaitsForActiveGroupsThenAcknowledgesEveryRank) {
  boost::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  alps::process_helper_mpi process(world, 1);
  alps::process_group active;
  EXPECT_FALSE(process.is_halting());
  if (world.rank() == 0) {
    active = process.allocate();
    process.halt();
    EXPECT_TRUE(process.is_halting());
  }
  world.barrier();
  // The scheduler cannot send shutdown while a worker group is active.
  EXPECT_FALSE(process.check_halted());
  world.barrier();
  if (world.rank() == 0)
    process.release(active);
  alps_test::expect_halted(process, world);
  EXPECT_TRUE(process.is_halting());
  EXPECT_TRUE(process.check_halted());
  process.halt();
  EXPECT_TRUE(process.check_halted());
}
