/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003 by Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id$ */

#ifndef FUNCTION_OBJECTS_HPP
#define FUNCTION_OBJECTS_HPP

#include <functional>

namespace boost
{
  // improved function objects taking optionally different argument types
  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct plus : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x + y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct minus : std::binary_function<Arg1, Arg2, Result>  {
    Result operator () (const Arg1& x, const Arg2& y) const { return x - y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct multiplies : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x * y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct divides : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x / y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct modulus : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x % y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=bool>
  struct logical_and : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x && y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=bool>
  struct logical_or : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x || y; }
  };

  // additional function objects for bit operations missing from the standard
  
  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct bit_and : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x & y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct bit_or : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x | y; }
  };

  template <class Arg1, class Arg2=Arg1, class Result=Arg1>
  struct bit_xor : std::binary_function<Arg1, Arg2, Result> {
    Result operator () (const Arg1& x, const Arg2& y) const { return x ^ y; }
  };

  template <class T>
  struct bit_not : std::unary_function<T, T> {
    T operator () (const T& x) const { return ~x; }
  };

} // namespace boost

#endif
