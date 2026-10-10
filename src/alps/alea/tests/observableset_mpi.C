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

#include <alps/alea.h>
#include <boost/mpi.hpp>
#include <boost/random.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <string>

TEST(ParallelObservableSet, PreservesSamplesAndDerivedStatisticsAcrossRanks) {
  boost::mpi::communicator world;
  ASSERT_GE(world.size(), 2);
  boost::mt19937 engine(2873u);
  boost::variate_generator<boost::mt19937 &, boost::uniform_real<>> random(engine,
                                                                           boost::uniform_real<>());
  alps::ObservableSet expected;
  expected << alps::RealObservable("observable a") << alps::RealObservable("observable b");
  for (int i = 0; i < (1 << 12); ++i) {
    expected["observable a"] << random();
    expected["observable b"] << random();
  }
  alps::RealObsevaluator a = expected["observable a"];
  alps::RealObsevaluator b = expected["observable b"];
  alps::RealObsevaluator ratio("ratio");
  ratio = a / b;
  expected.addObservable(ratio);

  // Relay the original seeded data through every rank, then back to the sender.
  // Every rank independently constructs the reference; it is not serialized.
  alps::ObservableSet received;
  if (world.rank() == 0) {
    world.send(1, 0, expected);
    world.recv(world.size() - 1, 0, received);
  } else {
    world.recv(world.rank() - 1, 0, received);
    world.send((world.rank() + 1) % world.size(), 0, received);
  }
  EXPECT_EQ(received.size(), expected.size());
  for (const std::string name : {"observable a", "observable b", "ratio"}) {
    SCOPED_TRACE("rank=" + std::to_string(world.rank()) + " observable=" + name);
    ASSERT_TRUE(received.has(name));
    const alps::RealObsevaluator actual = received[name];
    const alps::RealObsevaluator reference = expected[name];
    EXPECT_EQ(actual.count(), 4096u);
    EXPECT_DOUBLE_EQ(actual.mean(), reference.mean());
    EXPECT_DOUBLE_EQ(actual.error(), reference.error());
    EXPECT_EQ(actual.converged_errors(), reference.converged_errors());
    EXPECT_TRUE(std::isfinite(actual.mean()));
    EXPECT_GT(actual.error(), 0);
  }
}
