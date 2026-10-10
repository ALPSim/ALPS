/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2007 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "observable_checks.hpp"
#include <sstream>
#include <cmath>
#include <alps/alea.h>
#include <alps/parameter.h>
#include <boost/random.hpp>

TEST(AleaDetailedbinning, HistoricalNumericalScenarios)
{

  typedef boost::minstd_rand0 random_base_type;
  typedef boost::uniform_01<random_base_type> random_type;
  random_base_type random_int(1u); // Preserve the historical default seed.
  random_type random(random_int);

  alps::RealObservable obs_a("observable a");
  alps::RealObservable obs_b("observable b");
  alps::SimpleRealObservable obs_c("observable c");
  alps::RealObservable obs_d("observable d");

  // Preserve the former detailedbinning.input parameters.
  std::istringstream input("THERMALIZATION=500; STEPS=10000;");
  alps::Parameters parms(input);
  unsigned thermalization_steps=parms.value_or_default("THERMALIZATION",1000);
  unsigned number_of_steps=parms.value_or_default("STEPS",10000);

  for(unsigned i = 0; i < thermalization_steps; ++i){
    random();
    random();
  }

  for(unsigned i = 0; i < number_of_steps; ++i){
    obs_a << random();
    obs_b << random()+1;
    obs_c << 1.2;
    obs_d << 1.2;
  }

  EXPECT_NEAR(obs_a.error(0), 0.00290477, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 0, 10000u);
  EXPECT_NEAR(obs_a.error(1), 0.00290752, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 1, 5000u);
  EXPECT_NEAR(obs_a.error(2), 0.00287766, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 2, 2500u);
  EXPECT_NEAR(obs_a.error(3), 0.00286006, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 3, 1250u);
  EXPECT_NEAR(obs_a.error(4), 0.00277273, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 4, 625u);
  EXPECT_NEAR(obs_a.error(5), 0.00279057, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 5, 312u);
  EXPECT_NEAR(obs_a.error(6), 0.00282794, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 6, 156u);
  alps_test::expect_estimate(obs_a, {{0.499479, 0.00283, 5.00001e-07, 5.00001e-06}}, "obs_a");
  EXPECT_NEAR((obs_a).tau(), -0.0261, 5.00001e-05);
  EXPECT_NEAR(obs_b.error(0), 0.00288604, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 0, 10000u);
  EXPECT_NEAR(obs_b.error(1), 0.00287808, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 1, 5000u);
  EXPECT_NEAR(obs_b.error(2), 0.00287685, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 2, 2500u);
  EXPECT_NEAR(obs_b.error(3), 0.00286248, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 3, 1250u);
  EXPECT_NEAR(obs_b.error(4), 0.00275769, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 4, 625u);
  EXPECT_NEAR(obs_b.error(5), 0.00285328, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 5, 312u);
  EXPECT_NEAR(obs_b.error(6), 0.0027477, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 6, 156u);
  alps_test::expect_estimate(obs_b, {{1.49968, 0.00275, 5.00001e-06, 5.00001e-06}}, "obs_b");
  EXPECT_NEAR((obs_b).tau(), -0.0468, 5.00001e-05);
  alps_test::expect_estimate(obs_c, {{1.2, 0, 5.00001e-06, 0.}}, "obs_c");
  EXPECT_NEAR(obs_d.error(0), 0, 0.);
  EXPECT_EQ(obs_d.count() >> 0, 10000u);
  EXPECT_NEAR(obs_d.error(1), 0, 0.);
  EXPECT_EQ(obs_d.count() >> 1, 5000u);
  EXPECT_NEAR(obs_d.error(2), 0, 0.);
  EXPECT_EQ(obs_d.count() >> 2, 2500u);
  EXPECT_NEAR(obs_d.error(3), 0, 0.);
  EXPECT_EQ(obs_d.count() >> 3, 1250u);
  EXPECT_NEAR(obs_d.error(4), 0, 0.);
  EXPECT_EQ(obs_d.count() >> 4, 625u);
  EXPECT_NEAR(obs_d.error(5), 0, 0.);
  EXPECT_EQ(obs_d.count() >> 5, 312u);
  EXPECT_NEAR(obs_d.error(6), 0, 0.);
  EXPECT_EQ(obs_d.count() >> 6, 156u);
  alps_test::expect_estimate(obs_d, {{1.2, 0, 5.00001e-06, 0.}}, "obs_d");
  // Constant data has undefined normalized autocorrelation; the old renderer
  // substituted a display-only zero when the error was zero.
  EXPECT_TRUE(std::isnan(obs_d.tau()));

  alps::RealObsevaluator obseval_a(obs_a);
  alps::RealObsevaluator obseval_b(obs_b);
  alps::RealObsevaluator obseval_c;
  obseval_c = obseval_b / obseval_a;
  alps_test::expect_estimate(obseval_c, {{3.00282, 0.0157, 5.00001e-06, 5.00001e-05}}, "obseval_c");

}
