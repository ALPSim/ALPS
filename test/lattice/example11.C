/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2007 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/lattice.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main()
{

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    alps::Parameters p;
    std::cin >> p;
    p["LATTICE"] = "chain lattice";
    p["L"] = 4;
    alps::graph_helper<> g(p);
    std::vector<std::string> labels = g.momenta_labels();
    BOOST_FOREACH(std::string const& s, labels) std::cout << s << std::endl;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  exit(-1);
}
catch (...)
{
  std::cerr << "Caught unknown exception\n";
  exit(-2);
}
#endif
  return 0;
}
