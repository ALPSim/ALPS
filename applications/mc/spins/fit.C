/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <stdexcept>
#include "factory.h"
#include "fitting_scheduler.h"
#include "alps/osiris.h"

/**
 * Main function to be called to start the fitting process.
 * Starting from a parameter file, the scheduler is started.
 */
int main(int argc, char** argv)
{
  #ifndef BOOST_NO_EXCEPTIONS
    try {
  #endif
     // choose the right factory, depending on the input!
     return start_fitting(argc,argv,SpinFactory());
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
