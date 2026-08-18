/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: nobinning.h 3520 2009-12-11 16:49:53Z gamperl $ */


#ifndef ALPS_NUMERIC_SPECIAL_FUNCTIONS_HPP
#define ALPS_NUMERIC_SPECIAL_FUNCTIONS_HPP

#include <boost/accumulators/numeric/functional/vector.hpp>
#include <boost/math/special_functions.hpp>


namespace alps {
  namespace numeric {
    
    // define special powers
    template<class T> 
    inline T sq(T value) {
        using boost::numeric::operators::operator*;
        return value * value; 
    }

    template<class T>
    inline T cb(T value) { 
        using boost::numeric::operators::operator*;
        return value * value * value; 
    }

    template<class T>
    inline T cbrt(T value) { 
        return std::pow(value,1./3.); 
    }

    // define norm and r
    template <class T>
    inline T norm(T x, T y=T(), T z=T()) {
        using boost::numeric::operators::operator+;
        return (sq(x) + sq(y) + sq(z));
    }
    
    template <class T>
    inline T r(T x, T y=T(), T z=T()) {
        return std::sqrt(norm(x,y,z)); 
    }
  }
}


#endif
