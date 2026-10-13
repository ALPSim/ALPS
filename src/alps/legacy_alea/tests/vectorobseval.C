/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/alea.h>
#include <boost/random.hpp>
#include "observable_checks.hpp"
#include <iomanip>

TEST(AleaVectorobseval, HistoricalNumericalScenarios)
{

  typedef boost::minstd_rand0 random_base_type;
  typedef boost::uniform_01<random_base_type> random_type;
  random_base_type random_int(1u); // Preserve the historical default seed.
  random_type random(random_int);

  alps::RealVectorObservable obs_a("observable");

  for(int i=0; i < (1<<12); ++i) {
    std::valarray<double> obs(2);
    obs[0] = random();
    obs[1] = random() + 1;
    obs_a << obs;
  }
  EXPECT_NEAR(obs_a.error(0)[0], 0.00449726, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 0, 4096u);
  EXPECT_NEAR(obs_a.error(1)[0], 0.00447824, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 1, 2048u);
  EXPECT_NEAR(obs_a.error(2)[0], 0.00442771, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 2, 1024u);
  EXPECT_NEAR(obs_a.error(3)[0], 0.0043568, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 3, 512u);
  EXPECT_NEAR(obs_a.error(4)[0], 0.00432941, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 4, 256u);
  EXPECT_NEAR(obs_a.error(5)[0], 0.0045331, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 5, 128u);
  EXPECT_NEAR(obs_a.error(0)[1], 0.00453124, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 0, 4096u);
  EXPECT_NEAR(obs_a.error(1)[1], 0.00446591, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 1, 2048u);
  EXPECT_NEAR(obs_a.error(2)[1], 0.00441569, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 2, 1024u);
  EXPECT_NEAR(obs_a.error(3)[1], 0.00425621, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 3, 512u);
  EXPECT_NEAR(obs_a.error(4)[1], 0.00434305, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 4, 256u);
  EXPECT_NEAR(obs_a.error(5)[1], 0.00438853, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 5, 128u);
  alps_test::expect_estimate(obs_a, {{0.502618, 0.0045331, 5.00001e-07, 5.00001e-08}, {1.50074, 0.00438853, 5.00001e-06, 5.00001e-09}}, "obs_a");
  EXPECT_NEAR((obs_a).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((obs_a).tau()[1], -0.0309996, 5.00001e-08);

  alps::RealVectorObsevaluator veceval0(obs_a);
  alps::RealVectorObsevaluator veceval1(obs_a);
  std::valarray<double> vec(2);
  vec[0] = 3;
  vec[1] = -2;

  alps_test::expect_estimate(veceval0, {{0.502618, 0.00453255, 5.00001e-07, 5.00001e-09}, {1.50074, 0.00438799, 5.00001e-06, 5.00001e-09}}, "veceval0");
  EXPECT_NEAR((veceval0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0).tau()[1], -0.0309996, 5.00001e-08);

  alps_test::expect_estimate(veceval0 + 1.0, {{1.50262, 0.00453255, 5.00001e-06, 5.00001e-09}, {2.50074, 0.00438799, 5.00001e-06, 5.00001e-09}}, "veceval0 + 1.0");
  EXPECT_NEAR((veceval0 + 1.0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0 + 1.0).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(2.5 + veceval0, {{3.00262, 0.00453255, 5.00001e-06, 5.00001e-09}, {4.00074, 0.00438799, 5.00001e-06, 5.00001e-09}}, "2.5 + veceval0");
  EXPECT_NEAR((2.5 + veceval0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((2.5 + veceval0).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(veceval0 + veceval1, {{1.00524, 0.0090651, 5.00001e-06, 5.00001e-08}, {3.00149, 0.00877598, 5.00001e-06, 5.00001e-09}}, "veceval0 + veceval1");
  alps_test::expect_estimate(vec + veceval0, {{3.50262, 0.00453255, 5.00001e-06, 5.00001e-09}, {-0.499255, 0.00438799, 5.00001e-07, 5.00001e-09}}, "vec + veceval0");
  EXPECT_NEAR((vec + veceval0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((vec + veceval0).tau()[1], -0.0309996, 5.00001e-08);

  alps_test::expect_estimate(veceval0 - 1.0, {{-0.497382, 0.00453255, 5.00001e-07, 5.00001e-09}, {0.500745, 0.00438799, 5.00001e-07, 5.00001e-09}}, "veceval0 - 1.0");
  EXPECT_NEAR((veceval0 - 1.0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0 - 1.0).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(2.5 - veceval0, {{1.99738, 0.00453255, 5.00001e-06, 5.00001e-09}, {0.999255, 0.00438799, 5.00001e-07, 5.00001e-09}}, "2.5 - veceval0");
  EXPECT_NEAR((2.5 - veceval0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((2.5 - veceval0).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(veceval0 - veceval1, {{0, 0, 0., 0.}, {0, 0, 0., 0.}}, "veceval0 - veceval1");
  alps_test::expect_estimate(veceval0 - vec, {{-2.49738, 0.00453255, 5.00001e-06, 5.00001e-09}, {3.50074, 0.00438799, 5.00001e-06, 5.00001e-09}}, "veceval0 - vec");
  EXPECT_NEAR((veceval0 - vec).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0 - vec).tau()[1], -0.0309996, 5.00001e-08);

  alps_test::expect_estimate(veceval0 * 3.0, {{1.50785, 0.0135976, 5.00001e-06, 5.00001e-08}, {4.50223, 0.013164, 5.00001e-06, 5.00001e-07}}, "veceval0 * 3.0");
  EXPECT_NEAR((veceval0 * 3.0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0 * 3.0).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(1.3 * veceval0, {{0.653403, 0.00589231, 5.00001e-07, 5.00001e-09}, {1.95097, 0.00570439, 5.00001e-06, 5.00001e-09}}, "1.3 * veceval0");
  EXPECT_NEAR((1.3 * veceval0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((1.3 * veceval0).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(veceval0 * veceval1, {{0.252604, 0.00455558, 5.00001e-07, 5.00001e-09}, {2.25222, 0.0131702, 5.00001e-06, 5.00001e-08}}, "veceval0 * veceval1");
  alps_test::expect_estimate(veceval0 * vec, {{1.50785, 0.0135976, 5.00001e-06, 5.00001e-08}, {-3.00149, 0.00877598, 5.00001e-06, 5.00001e-09}}, "veceval0 * vec");
  EXPECT_NEAR((veceval0 * vec).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0 * vec).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(vec * veceval0, {{1.50785, 0.0135976, 5.00001e-06, 5.00001e-08}, {-3.00149, 0.00877598, 5.00001e-06, 5.00001e-09}}, "vec * veceval0");
  EXPECT_NEAR((vec * veceval0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((vec * veceval0).tau()[1], -0.0309996, 5.00001e-08);

  alps_test::expect_estimate(veceval0 / 3.0, {{0.167539, 0.00151085, 5.00001e-07, 5.00001e-09}, {0.500248, 0.00146266, 5.00001e-07, 5.00001e-09}}, "veceval0 / 3.0");
  EXPECT_NEAR((veceval0 / 3.0).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0 / 3.0).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(1.3 / veceval0, {{2.58646, 0.0233244, 5.00001e-06, 5.00001e-08}, {0.866237, 0.00253277, 5.00001e-07, 5.00001e-09}}, "1.3 / veceval0");
  alps_test::expect_estimate(veceval0 / veceval1, {{1, 0, 5.00001e-06, 0.}, {1, 0, 5.00001e-06, 0.}}, "veceval0 / veceval1");
  alps_test::expect_estimate(veceval0 / vec, {{0.167539, 0.00151085, 5.00001e-07, 5.00001e-09}, {-0.750372, 0.002194, 5.00001e-07, 5.00001e-07}}, "veceval0 / vec");
  EXPECT_NEAR((veceval0 / vec).tau()[0], 0.00800212, 5.00001e-09);
  EXPECT_NEAR((veceval0 / vec).tau()[1], -0.0309996, 5.00001e-08);
  alps_test::expect_estimate(vec / veceval0, {{5.96875, 0.0538255, 5.00001e-06, 5.00001e-08}, {-1.33267, 0.00389657, 5.00001e-06, 5.00001e-09}}, "vec / veceval0");

  alps_test::expect_estimate(pow(veceval0, 2), {{0.252604, 0.00455558, 5.00001e-07, 5.00001e-09}, {2.25222, 0.0131702, 5.00001e-06, 5.00001e-08}}, "pow(veceval0, 2)");

  alps::RealObsevaluator eval0(obs_a.slice(0));
  alps::RealObsevaluator eval1(obs_a.slice(1));
  alps_test::expect_estimate(eval0, {{0.502618, 0.00453, 5.00001e-07, 5.00001e-06}}, "eval0");
  EXPECT_NEAR((eval0).tau(), 0.008, 0.000500001);
  alps_test::expect_estimate(eval1, {{1.50074, 0.00439, 5.00001e-06, 5.00001e-06}}, "eval1");
  EXPECT_NEAR((eval1).tau(), -0.031, 0.000500001);
  alps_test::expect_estimate(eval0 * eval1, {{0.7543, 0.00727, 5.00001e-07, 5.00001e-06}}, "eval0 * eval1");

}
