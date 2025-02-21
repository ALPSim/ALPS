/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2008 by Synge Todo <wistaria@comp-phys.org>
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

/* $Id: parameters.C 2853 2008-06-17 13:59:59Z wistaria $ */

#include <alps/parameter/parameters.h>
#include <boost/config.hpp>
#include <boost/mpi.hpp>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

int main(int argc, char* argv[])
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  if (world.size() >= 2) {
    if (world.rank() == 0) {
      alps::Parameters params(std::cin);
      world.send(1, 0, params);
    } else if (world.rank() == 1) {
      alps::Parameters params;
      world.recv(0, 0, params);
      std::cout << params;
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e) {
  std::cerr << "Caught exception: " << e.what() << "\n";
  exit(-1);
}
catch (...) {
  std::cerr << "Caught unknown exception\n";
  exit(-2);
}
#endif
}
