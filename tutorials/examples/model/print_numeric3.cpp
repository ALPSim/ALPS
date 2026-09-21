/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*               2010-2010 by Ryo IGARASHI <rigarash@hosi.phys.s.u-tokyo.ac.jp>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <alps/model/hamiltonian_matrix.hpp>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  using namespace boost::numeric::ublas ;

  alps::Parameters parms;
  std::cin >> parms;

  alps::hamiltonian_matrix<compressed_matrix<double, row_major> > matrix(parms);

  if (matrix.uses_translation_invariance()) {
      std::cout << matrix.bloch_states_vector() << "\n";
  } else {
      std::cout << matrix.states_vector() << "\n";
  }
  std::cout << matrix.matrix() << "\n";

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
