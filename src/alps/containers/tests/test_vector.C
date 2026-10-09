/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/fixed_capacity_vector.h>
#include <cstddef>
#include <iostream>
#include <vector>

#include "test_main.h"

const std::size_t m = 10000; // number of trials
const std::size_t n = 16;    // max elements

int main() {
  test_main<alps::fixed_capacity_vector<double, n>,
            alps::fixed_capacity_vector<non_pod<double>, n>,
            std::vector<double> >(m, n);
  std::cout << "Test done successfully.\n";
  return 0;
}
