/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2004-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>,
*                            Ian McCulloch <ianmcc@physik.rwth-aachen.de>
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

template <typename GraphType>
void IterateOverBonds(const GraphType& graph)
{
  // iterate over all bonds in the graph.

  // The iterator type is obtained from alps::graph_traits
  typedef typename alps::graph_traits<GraphType>::bond_iterator bond_iterator;

  // Each bond has a source and target site, labelled by the site_descriptor type
  // we encountered in example5.
  typedef typename alps::graph_traits<GraphType>::site_descriptor site_descriptor;
  BOOST_STATIC_ASSERT((boost::is_convertible<site_descriptor, int>::value));

  // Each bond as an associated integer 'bond_type', acessed through
  // the bond_type_t property tag.
  // (NOTE however that it is actually read from the xml file
  // as an unsigned int, so don't use negative numbers!)
  // This gives access to the 'type' attribute of the <EDGE> xml element.
  typename alps::property_map<alps::bond_type_t, GraphType ,int>::const_type
    bond_type(get_or_default(alps::bond_type_t(), graph, 0));

  // determine the total number of bonds in the graph
  std::cout << "The graph has " << num_bonds(graph) << " bonds.\n";

  // bonds(graph) returns a [begin, end) pair of iterators over all bonds
  bond_iterator bond_it, bond_end;
  for (boost::tie(bond_it, bond_end) = bonds(graph); bond_it != bond_end;
       ++bond_it) 
  {
    // determine the source and target sites of the bond.
    site_descriptor source = boost::source(*bond_it, graph);
    site_descriptor target = boost::target(*bond_it, graph);

    int type = bond_type[*bond_it];
    
    std::cout << "The bond between site " << source 
              << " and site " << target
              << " has type " << type << std::endl;
  }
}

int main()
{

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    // read parameters
    alps::Parameters parameters;
    std::cin >> parameters;
    // create a graph factory with default graph type
    alps::graph_helper<> lattice(parameters);

    IterateOverBonds(lattice.graph());

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
