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
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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
