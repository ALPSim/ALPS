/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_NUMERIC_IS_NONZERO_HPP
#define ALPS_NUMERIC_IS_NONZERO_HPP

#include <alps/numeric/is_zero.hpp>

namespace alps { namespace numeric {

//
// is_nonzero
//

/// \brief checks if a number is not zero
/// in case of a floating point number, absolute values less than
/// epsilon (1e-50 by default) count as zero
/// \return returns true if the value is not zero

template<unsigned int N, class T>
inline bool is_nonzero(T x)
{ 
  return !is_zero<N>(x); 
}


template<class T>
inline bool is_nonzero(T x)
{ 
  return !is_zero(x); 
}


} } // end namespace alps::alea

#endif // ALPS_NUMERIC_IS_NONZERO_HPP
