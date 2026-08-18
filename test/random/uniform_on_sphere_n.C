/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/random/uniform_on_sphere_n.h>
#include <alps/utility/vectorio.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <vector>

int main() {
  static const int dim = 3;
  boost::mt19937 eng;
  alps::uniform_on_sphere_n<dim, double, std::vector<double> > dist;
  for (int i = 0; i < 100; ++i) {
    std::vector<double> r = dist(eng);
    std::cout << '(' << alps::write_vector(r, ",", 3) << ')' << std::endl;
  }
  return 0;
}
