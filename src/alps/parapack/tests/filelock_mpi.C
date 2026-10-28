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
#include <alps/parapack/mpi.h>
#include <gtest/gtest.h>
#include <memory>
#include <chrono>
#include <thread>
#include <stdexcept>
#include <string>

TEST(ParallelFileLock, ExcludesOtherRanksAndReleasesOnScopeExit) {
  alps::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  alps::testing::TemporaryDirectory directory;
  std::string filename = (directory.path() / "shared-checkpoint").string();
  alps::mpi::broadcast(world, filename, 0);
  const std::filesystem::path path(filename);
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

TEST(ParallelFileLock, BlockingLockRetriesUntilOwnerReleases) {
  alps::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  alps::testing::TemporaryDirectory directory;
  std::string filename = (directory.path() / "contended-checkpoint").string();
  alps::mpi::broadcast(world, filename, 0);
  alps::filelock lock{std::filesystem::path(filename)};
  if (world.rank() == 0)
    EXPECT_NO_THROW(lock.lock(0));
  world.barrier();
  if (world.rank() == 0) {
    int ready = 0, acquired = 0;
    alps::mpi::recv(world, 1, 0, ready);
    MPI_Request completed;
    MPI_Irecv(&acquired, 1, MPI_INT, 1, 1, world, &completed);
    // Keep the lock through at least one one-second retry. CTest bounds the
    // entire topology if a regression makes the blocking call hang.
    std::this_thread::sleep_for(std::chrono::seconds(2));
    int done = 0;
    MPI_Test(&completed, &done, MPI_STATUS_IGNORE);
    EXPECT_FALSE(done) << "contender returned before release";
    EXPECT_NO_THROW(lock.release());
    if (!done) MPI_Wait(&completed, MPI_STATUS_IGNORE);
    EXPECT_EQ(acquired, 1);
  } else if (world.rank() == 1) {
    alps::mpi::send(world, 0, 0, 1);
    EXPECT_NO_THROW(lock.lock());
    const int acquired = lock.locking() ? 1 : 0;
    alps::mpi::send(world, 0, 1, acquired);
    if (lock.locking())
      EXPECT_NO_THROW(lock.release());
  }
  world.barrier();
  EXPECT_FALSE(lock.locked());
  world.barrier();
}
