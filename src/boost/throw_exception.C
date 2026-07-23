/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <boost/throw_exception.hpp>

#ifdef BOOST_NO_EXCEPTIONS

#include <cstdlib>
#include <iostream>
#include <stdexcept>

// Implementation of boost::throw_exception: if exception handling is
// disabled, just print a message to std::cerr and abort.

void boost::throw_exception(std::exception const& e)
{
  std::cerr << "Exception occured: " << e.what() << "\n"; 
  std::abort();
}

#endif
