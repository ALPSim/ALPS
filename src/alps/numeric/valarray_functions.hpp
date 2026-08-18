/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: valarray_functions.h 3520 2009-12-11 16:49:53Z tamama $ */

#ifndef ALPS_NUMERIC_VALARRAY_FUNCTIONS_HPP
#define ALPS_NUMERIC_VALARRAY_FUNCTIONS_HPP



namespace alps { 
  namespace numeric {

    template <class T>
    std::ostream& operator<< (std::ostream &out, std::valarray<T> const & val)
    {
      std::copy(&const_cast<std::valarray<T>&>(val)[0],&const_cast<std::valarray<T>&>(val)[0]+val.size(),std::ostream_iterator<T>(out,"\t"));
      return out;
    }


  }
}

#endif // ALPS_NUMERIC_VALARRAY_FUNCTIONS_HPP




