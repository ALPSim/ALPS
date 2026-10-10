/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2013 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/random/random_choice.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

static const unsigned int n = 9;
static const unsigned int samples = 100000;

TEST(RandomChoice, SeededWeightsAndFrequencies) {
  std::cout << "number of bins = " << n << std::endl;
  std::cout << "number of samples = " << samples << std::endl;

  std::vector<double> weights(n);

  // random number generator
  typedef boost::mt19937 engine_type;
  typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<> > generator_type;
  engine_type eng(29411);
  generator_type rng(eng, boost::uniform_real<>());

  // generate weights
  std::generate(weights.begin(), weights.end(), rng);
  double tw = 0;
  for (unsigned int i = 0; i < n; ++i) tw += weights[i];

  // double-base version
  {
    alps::random_choice<generator_type> dist(weights);
    EXPECT_TRUE(dist.check(weights));

    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[dist(rng)];

    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / tw) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / tw) << "    \t"
                << (accum[i] / samples) << "    \t" << diff << "    \t"
                << sigma << "    \t" << (diff / sigma) << std::endl;
      EXPECT_LT(diff, 3 * sigma) << "seed=29411, bin=" << i;
    }
  }

  // integer-base version
  {
    alps::random_choice<engine_type> dist(weights);
    EXPECT_TRUE(dist.check(weights));

    std::vector<double> accum(n, 0);
    for (unsigned int t = 0; t < samples; ++t) ++accum[dist(eng)];

    std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
    for (unsigned int i = 0; i < n; ++i) {
      double diff = std::abs((weights[i] / tw) - (accum[i] / samples));
      double sigma = std::sqrt(accum[i]) / samples;
      std::cout << i << "\t" << (weights[i] / tw) << "    \t"
                << (accum[i] / samples) << "    \t" << diff << "    \t"
                << sigma << "    \t" << (diff / sigma) << std::endl;
      EXPECT_LT(diff, 3 * sigma) << "seed=29411, bin=" << i;
    }
  }
}
