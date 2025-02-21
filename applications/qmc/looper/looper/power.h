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

#ifndef LOOPER_POWER_H
#define LOOPER_POWER_H

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>

namespace looper {

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

//
// function power2, power3, power4
//

#ifndef BOOST_NO_SFINAE

template<typename T>
T power2(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(std::complex<T> const& t)
{ return power2(real(t)) + power2(imag(t)); }

template<typename T>
T power3(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power3(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power4(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(std::complex<T> const& t)
{ return power2(power2(t)); }

#else

template<typename T>
T power2(T const& t) { return t * t; }

template<typename T>
T power3(T const& t) { return t * t * t; }

template<typename T>
T power4(T const& t) { return power2(power2(t)); }

#endif

} // end namespace looper

#endif // LOOPER_POWER_H
