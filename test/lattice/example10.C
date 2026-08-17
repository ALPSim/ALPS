/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/lattice.h>
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

    // read parameters
    alps::Parameters parameters;
    std::cin >> parameters;
    // create a graph factory with default graph type
    alps::graph_helper<> h(parameters);
    // iterate over all cells and test the index function
    typedef alps::graph_helper<>::lattice_type lattice_type;
    lattice_type const& l(h.lattice());
    typedef alps::lattice_traits<lattice_type>::cell_iterator cell_iterator;
    cell_iterator it,last;
    int i=0;
    for (boost::tie(it,last) = alps::cells(l); it != last ; ++it, ++i)
      std::cout << "Cell " << i 
                << " has offset " << alps::write_vector(alps::offset(*it,l))
                << " and index " << alps::index(*it,l) << "\n";

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
