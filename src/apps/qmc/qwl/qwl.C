/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004 by Stefan Wessel <wessel@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "qwl_sse.h"
#include "alps/osiris/comm.h"

typedef alps::scheduler::SimpleMCFactory<QWL_SSE_Simulation> QWL_SSE_Factory;

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
   return alps::scheduler::start(argc,argv,QWL_SSE_Factory());
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
      comm_exit(true);
      return -1;
    }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
}

