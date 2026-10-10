// SPDX-License-Identifier: MIT
#pragma once

#include <alps/parapack/process.h>
#include <gtest/gtest.h>
#include <chrono>
#include <thread>

namespace alps_test {
inline void expect_halted(alps::process_helper_mpi &process,
                          const boost::mpi::communicator &world) {
  process.halt();
  const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
  bool halted = false;
  do {
    halted = process.check_halted();
    if (!halted)
      std::this_thread::yield();
  } while (!halted && std::chrono::steady_clock::now() < deadline);
  EXPECT_TRUE(halted) << "rank=" << world.rank();
  world.barrier();
}
} // namespace alps_test
