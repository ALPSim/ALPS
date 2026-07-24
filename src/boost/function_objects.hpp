/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003 by Matthias Troyer <troyer@comp-phys.org>
*
* SPDX-License-Identifier: MIT
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
