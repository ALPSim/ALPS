/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                       Matthias Troyer <troyer@itp.phys.ethz.ch>,
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


#ifndef ALPS_NUMERIC_ACCUMULATE_IF_HPP
#define ALPS_NUMERIC_ACCUMULATE_IF_HPP

#include <iostream>
#include <iterator>


namespace alps {
namespace numeric {


template <class InputIterator, class BinaryOperation, class BinaryPredicate, class T>
T
  accumulate_if
    (  InputIterator    first
    ,  InputIterator    last
    ,  T                init
    ,  BinaryOperation  binary_op
    ,  BinaryPredicate  pred
    )
  {
    while (first != last)
    {
      if (pred(*first))
        init = binary_op(init, *first);
      ++first;
    }
    return init;
  }


template <class InputIterator, class BinaryPredicate, class T>
T
  accumulate_if
    (  InputIterator    first
    ,  InputIterator    last
    ,  T                init
    ,  BinaryPredicate  pred
    )
  {
    return accumulate_if(first, last, init, std::plus<T>(), pred);
  } 


} // ending namespace numeric
} // ending namespace alps


#endif
