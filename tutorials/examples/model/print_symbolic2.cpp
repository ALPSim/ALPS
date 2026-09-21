/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/expression.h>

#include <alps/model/hamiltonian_matrix.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  using namespace boost::numeric::ublas ;

  alps::Parameters parms;
  std::cin >> parms;

  alps::hamiltonian_matrix<compressed_matrix<alps::Expression,row_major> > ham(parms);

  std::cout << ham.matrix() << "\n";

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
}

