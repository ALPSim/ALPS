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

/* $Id: simpleobseval.C 3654 2010-01-06 23:47:46Z troyer $ */

#include <alps/alea.h>
#include <alps/alea/mcdata.hpp>
#include <boost/random.hpp>
#include "observable_checks.hpp"
#include <iomanip>

TEST(AleaMcdata, HistoricalNumericalScenarios)
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
  for(int i=0; i < 4; ++i) {
    obs_d << random();
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
  alps_test::expect_estimate(obs_d, {{0.664328, 0.163898, 5.00001e-07, 5.00001e-07}}, "obs_d");

  alps::alea::mcdata<double> obseval_a(obs_a);
  alps_test::expect_estimate(obseval_a, {{0.497211, 0.00484938, 5.00001e-07, 5.00001e-09}}, "obseval_a");
  EXPECT_NEAR((obseval_a).tau(), 0.0679088, 5.00001e-08);
  EXPECT_EQ(obseval_a.count(), 4096u);
  EXPECT_EQ(obseval_a.bin_size(), 32u);
  EXPECT_EQ(obseval_a.bin_number(), 128u);

  alps::alea::mcdata<double> obseval_b(obs_b);
  alps_test::expect_estimate(obseval_b, {{1.50455, 0.00470255, 5.00001e-06, 5.00001e-09}}, "obseval_b");
  EXPECT_NEAR((obseval_b).tau(), 0.0490359, 5.00001e-08);
  EXPECT_EQ(obseval_b.count(), 4096u);
  EXPECT_EQ(obseval_b.bin_size(), 32u);
  EXPECT_EQ(obseval_b.bin_number(), 128u);
  alps::alea::mcdata<double> obseval_c(obs_c);
  alps_test::expect_estimate(obseval_c, {{1.49876, 0.00255591, 5.00001e-06, 5.00001e-09}}, "obseval_c");
  EXPECT_NEAR((obseval_c).tau(), -0.018531, 5.00001e-07);
  EXPECT_EQ(obseval_c.count(), 12288u);
  EXPECT_EQ(obseval_c.bin_size(), 128u);
  EXPECT_EQ(obseval_c.bin_number(), 96u);
  alps::alea::mcdata<double> obseval_d(obs_d);
  alps_test::expect_estimate(obseval_d, {{0.664328, 0.163898, 5.00001e-07, 5.00001e-07}}, "obseval_d");
  EXPECT_EQ(obseval_d.count(), 4u);
  EXPECT_EQ(obseval_d.bin_size(), 1u);
  EXPECT_EQ(obseval_d.bin_number(), 4u);

  alps::alea::mcdata<double> obseval_0 = obseval_b / obseval_a;
  alps_test::expect_estimate(obseval_0, {{3.02568, 0.0309813, 5.00001e-06, 5.00001e-08}}, "obseval_0");

  alps::alea::mcdata<double> obseval_1(obseval_b / obseval_a);
  alps_test::expect_estimate(obseval_1, {{3.02568, 0.0309813, 5.00001e-06, 5.00001e-08}}, "obseval_1");

  alps::alea::mcdata<double> obseval_2(obseval_b / obseval_a);
  alps_test::expect_estimate(obseval_2, {{3.02568, 0.0309813, 5.00001e-06, 5.00001e-08}}, "obseval_2");

  alps::alea::mcdata<double> obseval_3;
  obseval_3 = (obseval_b / obseval_a);
  alps_test::expect_estimate(obseval_3, {{3.02568, 0.0309813, 5.00001e-06, 5.00001e-08}}, "obseval_3");

  alps::alea::mcdata<double> obseval_5;
  obseval_5 << obseval_b << obseval_c;
  alps_test::expect_estimate(obseval_5, {{1.50021, 0.00246925, 5.00001e-06, 5.00001e-09}}, "obseval_5");
  EXPECT_NEAR((obseval_5).tau(), -0.0016393, 5.00001e-08);
  EXPECT_EQ(obseval_5.count(), 16384u);
  EXPECT_EQ(obseval_5.bin_size(), 128u);
  EXPECT_EQ(obseval_5.bin_number(), 128u);

  alps::alea::mcdata<double> obseval_6;
  obseval_6 << obseval_b << obseval_d;
  alps_test::expect_estimate(obseval_6, {{1.50373, 0.00470068, 5.00001e-06, 5.00001e-09}}, "obseval_6");
  EXPECT_EQ(obseval_6.count(), 4100u);
  EXPECT_EQ(obseval_6.bin_size(), 32u);
  EXPECT_EQ(obseval_6.bin_number(), 128u);

  alps::alea::mcdata<double> obseval_h = 1.0 / obseval_a;
  alps_test::expect_estimate(obseval_h, {{2.01122, 0.0196157, 5.00001e-06, 5.00001e-08}}, "obseval_h");
  EXPECT_NEAR((obseval_h).tau(), 0.0679088, 5.00001e-08);

  alps::alea::mcdata<double> obseval_i = 1.0 / obseval_b;
  alps_test::expect_estimate(obseval_i, {{0.664652, 0.00207741, 5.00001e-07, 5.00001e-09}}, "obseval_i");
  EXPECT_NEAR((obseval_i).tau(), 0.0490359, 5.00001e-08);

  alps::alea::mcdata<double> obseval_j;
  obseval_j = pow(obseval_a, 3.3);
  alps_test::expect_estimate(obseval_j, {{0.0996385, 0.00320685, 5.00001e-08, 5.00001e-09}}, "obseval_j");
  EXPECT_NEAR((obseval_j).tau(), 0.0679088, 5.00001e-08);

  alps::alea::mcdata<double> obseval_k = obseval_6 / obseval_6;
  alps_test::expect_estimate(obseval_k, {{1, 0, 5.00001e-06, 0.}}, "obseval_k");

}
