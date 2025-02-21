/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#ifndef OSIRIS_STD_SET_H
#define OSIRIS_STD_SET_H

#include <alps/config.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/impl.h>

#include <set>

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

/// deserialize a std::set container
template <class T, class Compare, class Allocator>
inline alps::IDump& operator >> (alps::IDump& dump,
                                   std::set<T,Compare,Allocator>& x)
{
  return alps::detail::loadSetLikeContainer(dump,x);
}

/// serialize a std::set container
template <class T, class Compare, class Allocator>
inline alps::ODump& operator << (alps::ODump& dump,
                                   const std::set<T,Compare,Allocator>& x)
{
  alps::detail::saveContainer(dump,x);
  return dump;
}          

/// deserialize a std::multiset container
template <class T, class Compare, class Allocator>
inline alps::IDump& operator >> (alps::IDump& dump,
                                   std::multiset<T,Compare,Allocator>& x)
{
  return alps::detail::loadSetLikeContainer(dump,x);
}

                          
/// serialize a std::multiset container
template <class T, class Compare, class Allocator>
inline alps::ODump& operator << (alps::ODump& dump,
                                   const std::multiset<T,Compare,Allocator>& x)
{
  return alps::detail::saveContainer(dump,x);
}          

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_std_MAP_H
