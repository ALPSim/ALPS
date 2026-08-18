/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef OSIRIS_STD_STRING_H
#define OSIRIS_STD_STRING_H

#include <alps/osiris/dump.h>

#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <string>

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

//=======================================================================
// string templates
//-----------------------------------------------------------------------

template <class charT, class traits, class Allocator>
inline alps::IDump& operator >> (alps::IDump& dump,
  std::basic_string<charT,traits,Allocator>& s)
{
  uint32_t sz(dump);
  if(sz) {
    charT* t = new charT[sz+1];
    dump.read_string(sz+1,t);
    if(t[sz]!=charT(0))
      boost::throw_exception(std::runtime_error("string on dump not terminating with '\\0'"));
        s=t;
    delete[] t;
    if(s.length()!=sz)
      boost::throw_exception(std::runtime_error("string on dump has incorrect length"));
  } else {
    s="";
  }
  return dump;
}


inline alps::IDump& operator >> (alps::IDump& dump, std::string& s)
{
  dump.read_string(s);
  return dump;
}


template <class charT, class traits, class Allocator>
inline alps::ODump& operator << (alps::ODump& dump,
  const std::basic_string<charT,traits,Allocator>& s)
{
 dump << uint32_t(s.size());
 if(s.size())
   dump.write_string(s.size()+1,s.c_str());
 return dump;
}

inline alps::ODump& operator << (alps::ODump& dump, const std::string& s)
{
  dump.write_string(s);
  return dump;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_STD_STRING_H
