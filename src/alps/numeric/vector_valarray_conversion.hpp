/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Maximilian Poprawe <poprawem@ethz.ch>
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

#ifndef ALPS_VECTOR_VALARRAY_CONVERSION
#define ALPS_VECTOR_VALARRAY_CONVERSION

#include <alps/utility/data.hpp>
#include <vector>
#include <valarray>
#include <algorithm>


namespace alps {

  template <class T>
  struct vector2valarray_type {
    typedef T type;
  };

  template <class T, class A>
  struct vector2valarray_type<std::vector<T, A> > {
    typedef std::valarray<T> type;
  };


  template <class T>
  struct valarray2vector_type {
    typedef T type;
  };

  template <class T>
  struct valarray2vector_type<std::valarray<T> > {
    typedef std::vector<T> type;
  };


  namespace numeric {

    template <class T>
    T valarray2vector (const T& value) {return value;}

    template<class T>
    std::vector<T> valarray2vector(std::valarray<T> const & from)
    {
      std::vector<T> to;
      to.reserve(from.size());
      std::copy(data(from),data(from)+from.size(),std::back_inserter(to));
      return to;
    }



    template <class T>
    T vector2valarray (const T& value) {return value;}

    template<class T>
    std::valarray<T> vector2valarray(std::vector<T> const & from)
    {
      std::valarray<T> to(from.size());
      std::copy(from.begin(),from.end(),data(to));
      return to;
    }

    template<class T1, class T2>
    std::valarray<T2> vector2valarray(std::vector<T1> const & from)
    {
      std::valarray<T2> to(from.size());
      std::copy(from.begin(),from.end(),data(to));
      return to;
    }

  }
}

#endif
