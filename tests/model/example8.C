/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/parameter.h>
#include <iostream>
#include <string>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    alps::ParameterList parms;
    std::cin >> parms;
    for (int i=0;i<parms.size();++i) {
      alps::ModelLibrary models(parms[i]);
      alps::graph_helper<> lattice(parms[i]);
      alps::HamiltonianDescriptor<short> ham(models.get_hamiltonian(lattice,parms[i]));
      parms[i].copy_undefined(ham.default_parameters());
      ham.set_parameters(parms[i]);
      if (has_sign_problem(ham,lattice,parms[i]))
        std::cout << "Model " << i+1 << " has a sign problem.\n";
      else
        std::cout << "Model " << i+1 << " has no sign problem.\n";
    }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
