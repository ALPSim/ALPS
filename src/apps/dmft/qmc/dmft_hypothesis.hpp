// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT

#ifndef ALPS_DMFT_HYPOTHESIS_HPP
#define ALPS_DMFT_HYPOTHESIS_HPP

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <string>

namespace mc_hypothesis {
// Compares Monte Carlo averages with exact values in units of their error
// bars, as tests/cli/mc_hypothesis.py does (Wallerberger and Gull,
// Phys. Rev. E 96, 053303 (2017), Sec. II A): z = (mean - exact) / error is
// taken as standard normal, and the `size` comparisons of one family share
// the level alpha (Bonferroni), which holds however they are correlated.
class family {
public:
  explicit family(int size, double alpha = 1e-3) : size_(size), threshold_(alpha / size) {}
  ~family() { EXPECT_EQ(count_, size_) << "comparisons in the family"; }

  void expect(std::string const & name, double mean, double error, double exact) {
    ++count_;
    if (!(std::isfinite(error) && error > 0)) {
      ADD_FAILURE() << name << " has error " << error;
      return;
    }
    const double z = (mean - exact) / error, p = std::erfc(std::abs(z) / std::sqrt(2.));
    std::cout << name << " = " << mean << " +- " << error << ", exact " << exact << ", z = " << z << '\n';
    EXPECT_GE(p, threshold_) << name << " = " << mean << " +- " << error << " differs from exact " << exact
                             << " by " << z << " error bars";
  }

private:
  int size_, count_ = 0;
  double threshold_;
};
}

#endif
