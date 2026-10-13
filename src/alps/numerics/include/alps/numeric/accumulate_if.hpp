/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                       Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
