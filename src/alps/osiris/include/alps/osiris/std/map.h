/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef OSIRIS_std_MAP_H
#define OSIRIS_std_MAP_H

#include <alps/config.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/impl.h>
#include <alps/osiris/std/pair.h>

#include <map>

namespace alps {

/// deserialize a std::map container
template <class Key, class T, class Compare, class Allocator>
inline alps::IDump& operator >> (alps::IDump& dump,
                                   std::map<Key,T,Compare,Allocator>& x)
{
  x=std::map<Key,T,Compare,Allocator>();
  uint32_t n(dump);
  Key k;
  T v;
  while (n--)
    {
      dump >> k >> v;
      x[k]=v;
    }
  
  return dump;
}

/// serialize a std::map container
template <class Key, class T, class Compare, class Allocator>
inline alps::ODump& operator << (alps::ODump& dump,
                                   const std::map<Key,T,Compare,Allocator>& x)
{
  alps::detail::saveContainer(dump,x);
  return dump;
}          

/// deserialize a std::multimap container
template <class Key, class T, class Compare, class Allocator>
inline alps::IDump& operator >> (alps::IDump& dump,
                                   std::multimap<Key,T,Compare,Allocator>& x)
{
  return alps::detail::loadSetLikeContainer(dump,x);
}

                          
/// serialize a std::multimap container
template <class Key, class T, class Compare, class Allocator>
inline alps::ODump& operator << (alps::ODump& dump,
                                   const std::multimap<Key,T,Compare,Allocator>& x)
{
  return alps::detail::saveContainer(dump,x);
}          

} // end namespace alps

#endif // OSIRIS_std_MAP_H
