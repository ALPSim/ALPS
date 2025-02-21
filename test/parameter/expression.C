/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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
