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
   // The scheduler prints the help text. Detect the same query here, with the
   // rules of alps::handle_cli_information (cli.cpp): nothing after "--" counts.
   bool help = false;
   for (int i = 1; i < argc && std::strcmp(argv[i], "--") != 0; ++i)
     if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0)
       help = true;
   if (help) {
     alps::cli_mpi_guard mpi(argc, argv); // keeps the rank available after the scheduler's help
     const int result = alps::scheduler::start(argc,argv,SpinFactory(), "spinmc");
     if (alps::cli_is_master())
       std::cout << "\nCouplings: H = -sum s_i.J.s_j over bonds, so J>0 is ferromagnetic. Spins have\n"
                 << "length S (default 1). CONVENTION=\"quantum\" flips the sign (J>0 is\n"
                 << "antiferromagnetic), sets the length to sqrt(S(S+1)) with default S=1/2, and\n"
                 << "does not allow a magnetic field.\n";
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
