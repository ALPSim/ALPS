/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: example4.C 730 2004-04-02 18:20:12Z troyer $ */

#include <alps/lattice.h>
#include <alps/parameter_input.hpp>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main(int argc, char ** argv)
{

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    alps::Parameters parameters;
    switch(argc) {
       case 1:
         parameters = alps::read_parameter_text(std::cin);
         break;
       case 2:
         parameters = alps::read_parameter_file(argv[1]);
         break;
       default:
         std::cerr << "Usage: " << argv[0] << " [parameterfile]\n";
    } 
   
    // create a graph factory with default graph type
    alps::graph_helper<> lattice(parameters);
    // write the graph created from the input in XML
    std::cout << lattice.graph();

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
}
