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
  // Pin the historical observableset_mpi.op-2 values from ALPS at
  // 64d5c09e0, independently of the freshly computed round-trip reference.
  // Tolerances are half the final printed digit (plus roundoff).
  const auto& actual_a = received.get<alps::RealObservable>("observable a");
  const auto& actual_b = received.get<alps::RealObservable>("observable b");
  const alps::RealObsevaluator actual_ratio = received["ratio"];
  EXPECT_NEAR(actual_a.mean(), 0.5009, 0.0000500001);
  EXPECT_NEAR(actual_b.mean(), 0.504867, 0.000000500001);
  EXPECT_NEAR(actual_a.error(), 0.00482, 0.00000500001);
  EXPECT_NEAR(actual_b.error(), 0.00483, 0.00000500001);
  EXPECT_NEAR(actual_a.tau(), 0.0622, 0.0000500001);
  EXPECT_NEAR(actual_b.tau(), 0.0654, 0.0000500001);
  EXPECT_NEAR(actual_ratio.mean(), 0.992046, 0.000000500001);
  EXPECT_NEAR(actual_ratio.error(), 0.0139, 0.0000500001);
  const double errors_a[] = {0.00454298, 0.00450135, 0.00451505, 0.00452936, 0.00444092, 0.00481723};
  const double errors_b[] = {0.00454103, 0.00454096, 0.00448772, 0.00449454, 0.00460777, 0.00482904};
  for (unsigned level = 0; level < 6; ++level) {
    SCOPED_TRACE("binning level=" + std::to_string(level));
    EXPECT_NEAR(actual_a.error(level), errors_a[level], 0.00000000500001);
    EXPECT_NEAR(actual_b.error(level), errors_b[level], 0.00000000500001);
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
