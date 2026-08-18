/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "WRun.h"
#include <alps/osiris/comm.h>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
   return alps::scheduler::start(argc,argv,alps::scheduler::SimpleMCFactory<WRun>());
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
