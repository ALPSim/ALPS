/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_NUMERIC_UPDATE_MINMAX_HPP
#define ALPS_NUMERIC_UPDATE_MINMAX_HPP

#include <alps/type_traits/slice.hpp>

namespace alps { namespace numeric {

template <class T>
void update_max(T& lhs, T const& rhs)
{
  for (typename slice_index<T>::type it = slices(lhs).first; 
       it < slices(lhs).second && it < slices(rhs).second; ++it)
    if (slice_value(lhs,it) < slice_value(rhs,it))
      slice_value(lhs,it) = slice_value(rhs,it);
}

template <class T>
void update_min(T& lhs, T const& rhs)
{
  for (typename slice_index<T>::type it = slices(lhs).first; 
       it < slices(lhs).second && it < slices(rhs).second; ++it)
    if (slice_value(rhs,it) < slice_value(lhs,it))
      slice_value(lhs,it) = slice_value(rhs,it);
}



} } // end namespace alps::numeric

#endif // ALPS_NUMERIC_UPDATE_MINMAX_HPP
