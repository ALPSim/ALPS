/***************************************************************************
* PALM++ library
*
* boost/throw_exception.C  Implementation of boost::throw_exception
*
* $Id$
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>,
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
**************************************************************************/

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
