/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2008 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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

TEST(AleaSimpleobseval, HistoricalNumericalScenarios)
{

  typedef boost::minstd_rand0 random_base_type;
  typedef boost::uniform_01<random_base_type> random_type;
  random_base_type random_int(1u); // Preserve the historical default seed.
  random_type random(random_int);

  alps::RealObservable obs_a("observable a");
  alps::RealObservable obs_b("observable b");
  alps::RealObservable obs_c("observable c");
  alps::RealObservable obs_d("observable d");

  for(int i=0; i < (1<<12); ++i) {
    obs_a << random();
    obs_b << random()+1;
    obs_c << random()+1;
    obs_c << random()+1;
    obs_c << random()+1;
  }

  EXPECT_NEAR(obs_a.error(0), 0.00455021, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 0, 4096u);
  EXPECT_NEAR(obs_a.error(1), 0.00455922, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 1, 2048u);
  EXPECT_NEAR(obs_a.error(2), 0.00457997, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 2, 1024u);
  EXPECT_NEAR(obs_a.error(3), 0.00463017, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 3, 512u);
  EXPECT_NEAR(obs_a.error(4), 0.00479441, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 4, 256u);
  EXPECT_NEAR(obs_a.error(5), 0.00484938, 5.00001e-09);
  EXPECT_EQ(obs_a.count() >> 5, 128u);
  alps_test::expect_estimate(obs_a, {{0.497211, 0.00485, 5.00001e-07, 5.00001e-06}}, "obs_a");
  EXPECT_NEAR((obs_a).tau(), 0.0679, 5.00001e-05);
  EXPECT_NEAR(obs_b.error(0), 0.00448764, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 0, 4096u);
  EXPECT_NEAR(obs_b.error(1), 0.00457509, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 1, 2048u);
  EXPECT_NEAR(obs_b.error(2), 0.00454425, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 2, 1024u);
  EXPECT_NEAR(obs_b.error(3), 0.00455364, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 3, 512u);
  EXPECT_NEAR(obs_b.error(4), 0.00459102, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 4, 256u);
  EXPECT_NEAR(obs_b.error(5), 0.00470255, 5.00001e-09);
  EXPECT_EQ(obs_b.count() >> 5, 128u);
  alps_test::expect_estimate(obs_b, {{1.50455, 0.0047, 5.00001e-06, 5.00001e-05}}, "obs_b");
  EXPECT_NEAR((obs_b).tau(), 0.049, 0.000500001);
  EXPECT_NEAR(obs_c.error(0), 0.00260463, 5.00001e-09);
  EXPECT_EQ(obs_c.count() >> 0, 12288u);
  EXPECT_NEAR(obs_c.error(1), 0.00261035, 5.00001e-09);
  EXPECT_EQ(obs_c.count() >> 1, 6144u);
  EXPECT_NEAR(obs_c.error(2), 0.00263224, 5.00001e-09);
  EXPECT_EQ(obs_c.count() >> 2, 3072u);
  EXPECT_NEAR(obs_c.error(3), 0.00257175, 5.00001e-09);
  EXPECT_EQ(obs_c.count() >> 3, 1536u);
  EXPECT_NEAR(obs_c.error(4), 0.00260109, 5.00001e-09);
  EXPECT_EQ(obs_c.count() >> 4, 768u);
  EXPECT_NEAR(obs_c.error(5), 0.00263753, 5.00001e-09);
  EXPECT_EQ(obs_c.count() >> 5, 384u);
  EXPECT_NEAR(obs_c.error(6), 0.00255591, 5.00001e-09);
  EXPECT_EQ(obs_c.count() >> 6, 192u);
  alps_test::expect_estimate(obs_c, {{1.49876, 0.00256, 5.00001e-06, 5.00001e-06}}, "obs_c");
  EXPECT_NEAR((obs_c).tau(), -0.0185, 5.00001e-05);
  EXPECT_EQ(obs_d.count(), 0u);

  alps::RealObsevaluator obseval_a(obs_a);
  alps_test::expect_estimate(obseval_a, {{0.497211, 0.00485, 5.00001e-07, 5.00001e-06}}, "obseval_a");
  EXPECT_NEAR((obseval_a).tau(), 0.0679, 5.00001e-05);
  EXPECT_EQ(obseval_a.count(), 4096u);
  EXPECT_EQ(obseval_a.bin_size(), 32u);
  EXPECT_EQ(obseval_a.bin_number(), 128u);

  auto scaled = obseval_a;
  scaled *= 2u;
  EXPECT_DOUBLE_EQ(scaled.mean(), 2 * obseval_a.mean());
  EXPECT_DOUBLE_EQ(scaled.error(), 2 * obseval_a.error());
  scaled /= 2u;
  EXPECT_DOUBLE_EQ(scaled.mean(), obseval_a.mean());
  EXPECT_DOUBLE_EQ(scaled.error(), obseval_a.error());

  alps::RealObsevaluator obseval_b(obs_b);
  alps_test::expect_estimate(obseval_b, {{1.50455, 0.0047, 5.00001e-06, 5.00001e-05}}, "obseval_b");
  EXPECT_NEAR((obseval_b).tau(), 0.049, 0.000500001);
  EXPECT_EQ(obseval_b.count(), 4096u);
  EXPECT_EQ(obseval_b.bin_size(), 32u);
  EXPECT_EQ(obseval_b.bin_number(), 128u);
  alps::RealObsevaluator obseval_c(obs_c);
  alps_test::expect_estimate(obseval_c, {{1.49876, 0.0027, 5.00001e-06, 5.00001e-05}}, "obseval_c");
  EXPECT_NEAR((obseval_c).tau(), 0.0381, 5.00001e-05);
  EXPECT_EQ(obseval_c.count(), 12288u);
  EXPECT_EQ(obseval_c.bin_size(), 128u);
  EXPECT_EQ(obseval_c.bin_number(), 96u);
  alps::RealObsevaluator obseval_d(obs_d);
  EXPECT_EQ(obseval_d.count(), 0u);
  EXPECT_EQ(obseval_d.bin_size(), 0u);
  EXPECT_EQ(obseval_d.bin_number(), 0u);

  alps::RealObsevaluator obseval_0 = obseval_b / obseval_a;
  EXPECT_EQ(obseval_0.name(), "(observable b) / (observable a)");
  alps_test::expect_estimate(obseval_0, {{3.02568, 0.031, 5.00001e-06, 0.000500001}}, "obseval_0");

  alps::RealObsevaluator obseval_1(obseval_b / obseval_a);
  EXPECT_EQ(obseval_1.name(), "(observable b) / (observable a)");
  alps_test::expect_estimate(obseval_1, {{3.02568, 0.031, 5.00001e-06, 0.000500001}}, "obseval_1");

  alps::RealObsevaluator obseval_2(obseval_b / obseval_a, "obseval_e");
  EXPECT_EQ(obseval_2.name(), "obseval_e");
  alps_test::expect_estimate(obseval_2, {{3.02568, 0.031, 5.00001e-06, 0.000500001}}, "obseval_2");

  alps::RealObsevaluator obseval_3("obseval_f");
  obseval_3 = (obseval_b / obseval_a);
  EXPECT_EQ(obseval_3.name(), "obseval_f");
  alps_test::expect_estimate(obseval_3, {{3.02568, 0.031, 5.00001e-06, 0.000500001}}, "obseval_3");

  alps::RealObsevaluator obseval_4("obseval_g");
  obseval_4 << (obseval_b / obseval_a);
  EXPECT_EQ(obseval_4.name(), "obseval_g");
  alps_test::expect_estimate(obseval_4, {{3.02568, 0.031, 5.00001e-06, 0.000500001}}, "obseval_4");

  EXPECT_EQ(obseval_4.count(), 4096u);
  EXPECT_EQ(obseval_4.bin_size(), 32u);
  EXPECT_EQ(obseval_4.bin_number(), 128u);

  alps::RealObsevaluator obseval_5("obseval_5");
  obseval_5 << obseval_b << obseval_c;
  EXPECT_EQ(obseval_5.name(), "obseval_5");
  alps_test::expect_estimate(obseval_5, {{1.50021, 0.00247, 5.00001e-06, 5.00001e-06}}, "obseval_5");
  EXPECT_NEAR((obseval_5).tau(), 0.101, 0.000500001);
  EXPECT_EQ(obseval_5.count(), 16384u);
  EXPECT_EQ(obseval_5.bin_size(), 128u);
  EXPECT_EQ(obseval_5.bin_number(), 128u);

  alps::RealObsevaluator obseval_6("obseval_6");
  obseval_6 << obseval_b << obseval_d;
  EXPECT_EQ(obseval_6.name(), "obseval_6");
  alps_test::expect_estimate(obseval_6, {{1.50455, 0.0047, 5.00001e-06, 5.00001e-05}}, "obseval_6");
  EXPECT_NEAR((obseval_6).tau(), 0.049, 0.000500001);
  EXPECT_EQ(obseval_6.count(), 4096u);
  EXPECT_EQ(obseval_6.bin_size(), 32u);
  EXPECT_EQ(obseval_6.bin_number(), 128u);

  alps::RealObsevaluator obseval_h = 1.0 / obseval_a;
  alps_test::expect_estimate(obseval_h, {{2.01122, 0.0196, 5.00001e-06, 5.00001e-05}}, "obseval_h");

  alps::RealObsevaluator obseval_i = 1.0 / obseval_b;
  alps_test::expect_estimate(obseval_i, {{0.664652, 0.00208, 5.00001e-07, 5.00001e-06}}, "obseval_i");

  alps::RealObsevaluator obseval_j("pow(observable a,3.3)");
  obseval_j = pow(obseval_a, 3.3);
  EXPECT_EQ(obseval_j.name(), "pow(observable a,3.3)");
  alps_test::expect_estimate(obseval_j, {{0.0996385, 0.00321, 5.00001e-08, 5.00001e-06}}, "obseval_j");

  alps::RealObsevaluator obseval_k = obseval_6 / obseval_6;
  alps_test::expect_estimate(obseval_k, {{1, 0, 5.00001e-06, 0.}}, "obseval_k");

}
