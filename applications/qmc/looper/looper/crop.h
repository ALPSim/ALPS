/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_CROP_H
#define LOOPER_CROP_H

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>

namespace looper {

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

//
// function crop_0, crop_01
//

template<typename T>
T crop_0(T x, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_0(T const& x, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_01(T x, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

template<typename T>
T crop_01(T const& x, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

} // end namespace looper

#endif // LOOPER_CROP_H
