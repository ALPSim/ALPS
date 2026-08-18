/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "factory.h"
#include <alps/scheduler.h>
#include <complex>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

   return alps::scheduler::start(argc,argv,FullDiagFactory());

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
