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

#include <alps/parapack/filelock.h>
#include <alps/testing/temporary_directory.hpp>
#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <string>

TEST(ParallelFileLock, ExcludesOtherRanksAndReleasesOnScopeExit) {
  boost::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  alps::testing::TemporaryDirectory directory;
  std::string filename = (directory.path() / "shared-checkpoint").string();
  boost::mpi::broadcast(world, filename, 0);
  const boost::filesystem::path path(filename);
  for (int owner = 0; owner < world.size(); ++owner) {
    SCOPED_TRACE("owner=" + std::to_string(owner));
    auto lock = std::make_unique<alps::filelock>(path);
    EXPECT_FALSE(lock->locking());
    EXPECT_FALSE(lock->locked());
    world.barrier();
    if (world.rank() == owner) {
      EXPECT_NO_THROW(lock->lock(0));
      EXPECT_TRUE(lock->locking());
    }
    world.barrier();
    EXPECT_TRUE(lock->locked());
    if (world.rank() != owner) {
      // One attempt only: a regression must fail rather than block indefinitely.
      EXPECT_THROW(lock->lock(0), std::logic_error);
      EXPECT_FALSE(lock->locking());
    }
    world.barrier();
    // The next rank can acquire after either explicit or RAII release.
    if (world.rank() == owner && owner % 2 == 0) {
      EXPECT_NO_THROW(lock->release());
      EXPECT_FALSE(lock->locking());
    }
    lock.reset();
    world.barrier();
    alps::filelock unlocked(path);
    EXPECT_FALSE(unlocked.locked());
    world.barrier();
  }
  // Keep the root-owned shared directory alive until all ranks have finished.
  world.barrier();
}
