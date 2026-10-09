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

/* $Id$ */

#include <alps/expression.h>

#include <boost/throw_exception.hpp>
#include <iostream>
#include <stdexcept>

void print_simplified(const std::string& e)
{
    alps::Expression x(e);
    x.simplify();
    std::cout << x << "\n";
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    print_simplified("3*(a*b)*2+5*(x*y)^2*x*3");
    print_simplified("3*a*b+5*a-3*a*b+2.5*a");
    print_simplified("3*a*b+5*a+3*a*b+2.5*a");

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
