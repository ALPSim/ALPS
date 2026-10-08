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
#include <iostream>

namespace mpi = boost::mpi;

#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>

TEST(ParallelVectors, UnevenDistributionAndCollection) {
  mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  boost::mt19937 engine(5489);
  std::vector<int> counts, offsets;
  int n = 0;
  for (int rank = 0; rank < world.size(); ++rank) {
    counts.push_back((engine() & 3) + 1);
    offsets.push_back(n);
    n += counts.back();
  }
  std::vector<uint32_t> expected;
  for (int i = 0; i < n; ++i) expected.push_back(engine());
  std::vector<uint32_t> source = world.rank() == 0 ? expected : std::vector<uint32_t>{};
  std::vector<uint32_t> distributed;
  alps::distribute_vector(world, counts, offsets, source, distributed);
  const auto first = expected.begin() + offsets[world.rank()];
  EXPECT_EQ(distributed, std::vector<uint32_t>(first, first + counts[world.rank()]))
      << "rank=" << world.rank();
  std::vector<uint32_t> collected(world.rank() == 0 ? n : 0);
  alps::collect_vector(world, counts, offsets, distributed, collected);
  if (world.rank() == 0) EXPECT_EQ(collected, expected);
}
