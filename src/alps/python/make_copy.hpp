/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Matthias Troyer <troyer@comp-phys.org>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_PYTHON_MAKE_COPY_HPP
#define ALPS_PYTHON_MAKE_COPY_HPP

#include <boost/python/dict.hpp>
namespace alps { namespace python {

template<class T>
T make_copy(T const& x, boost::python::dict const& ) { return x; } 
 
} } // end namespace alps::python

#endif // ALPS_PYTHON_MAKE_COPY_HPP
