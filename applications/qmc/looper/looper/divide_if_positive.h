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

#ifndef LOOPER_DIVIDE_IF_POSITIVE_H
#define LOOPER_DIVIDE_IF_POSITIVE_H

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>

namespace looper {

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

//
// dip (divide_if_positive)
//

template<typename T, typename U>
T dip(T x, U y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T x, U const& y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T const& x, U y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T const& x, U const& y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T divide_if_positive(T x, U y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

template<typename T, typename U>
T divide_if_positive(T x, U const& y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

template<typename T, typename U>
T divide_if_positive(T const& x, U y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

template<typename T, typename U>
T divide_if_positive(T const& x, U const& y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

} // end namespace looper

#endif // LOOPER_DIVIDE_IF_POSITIVE_H
