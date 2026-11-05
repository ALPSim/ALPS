/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: parameters.C 2853 2008-06-17 13:59:59Z wistaria $ */

#include <alps/parameter/parameters.h>
#include <alps/parapack/mpi.h>
#include <gtest/gtest.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

TEST(ParallelLegacyParameters, PreservesKeysValuesExpressionsAndOrder) {
  alps::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  std::istringstream input("L=10; M=1; T=0.1; beta=\"1/T\";");
  const alps::Parameters original(input);
  alps::Parameters restored;
  if (world.rank() == 0) {
    alps::mpi::send(world, 1, 0, original);
    alps::mpi::recv(world, world.size() - 1, 0, restored);
  } else {
    alps::mpi::recv(world, world.rank() - 1, 0, restored);
    alps::mpi::send(world, (world.rank() + 1) % world.size(), 0, restored);
  }
  const std::vector<std::pair<std::string, std::string>> expected{
      {"L", "10"}, {"M", "1"}, {"T", "0.1"}, {"beta", "1/T"}};
  ASSERT_EQ(restored.size(), expected.size());
  auto actual = restored.begin();
  for (const auto &entry : expected) {
    SCOPED_TRACE("rank=" + std::to_string(world.rank()) + " key=" + entry.first);
    EXPECT_EQ(actual->key(), entry.first);
    EXPECT_EQ(std::string(actual->value().c_str()), entry.second);
    ++actual;
  }
}
