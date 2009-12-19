/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
