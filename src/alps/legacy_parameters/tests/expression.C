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

#include <alps/expression.h>

#include <boost/throw_exception.hpp>
#include <iostream>
#include <stdexcept>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  alps::Parameters parms;
  std::string str;
  while (std::getline(std::cin, str) && str.size() && str[0] != '%')
    parms.push_back(alps::Parameter(str));
  std::cout << "Parameters:\n" << parms << std::endl;
  
  alps::ParameterEvaluator eval(parms);
  while (std::cin) {
    alps::Expression expr(std::cin);
    if (!expr.can_evaluate(eval))
      std::cout << "Cannot evaluate [" << expr << "]." << std::endl;
    else 
      std::cout << "The value of [" << expr << "] is "
                << alps::evaluate<double>(expr, eval) << std::endl;
    char c;
    std::cin >> c;
    if (c!=',')
      break;
  }

  while (std::cin) {
    std::string v;
    std::cin >> v;
    if (v.empty()) break;
    if (!alps::can_evaluate(v, parms))
      std::cout << "Cannot evaluate [" << v << "]." << std::endl;
    else 
      std::cout << "The value of [" << v << "] is "
                << alps::evaluate<double>(v, parms) << std::endl;
  }

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
