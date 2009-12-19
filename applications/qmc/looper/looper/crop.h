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
