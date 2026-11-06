/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2003 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/cli.hpp>
#include <alps/scheduler.h>
#include <alps/osiris/comm.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include "factory.h"

/**
 * Main function to be called to start the fitting process.
 * Starting from a parameter file, the scheduler is started.
 */

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
   const bool help = std::any_of(argv + 1, argv + argc, [](const char* argument) {
     return std::strcmp(argument, "--help") == 0 || std::strcmp(argument, "-h") == 0;
   });
   if (help) {
     alps::cli_mpi_guard mpi(argc, argv); // keeps the rank available after the scheduler's help
     const int result = alps::scheduler::start(argc,argv,SpinFactory(), "spinmc");
     if (alps::cli_is_master())
       std::cout << "\nJ>0 is ferromagnetic, H = -J sum s_i.s_j. CONVENTION=\"quantum\" makes it\n"
                 << "antiferromagnetic and sets the spin length to sqrt(S(S+1)).\n";
     return result;
   }
   return alps::scheduler::start(argc,argv,SpinFactory(), "spinmc");
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
      alps::comm_exit(true);
      return -1;
    }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
}
