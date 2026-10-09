/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2011 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef OSIRIS_BOOST_UBLAS_H
#define OSIRIS_BOOST_UBLAS_H

// #include <palm/config.h>
#include <boost/numeric/ublas/vector.hpp>
#include <alps/osiris/std/impl.h>

/// deserialize a boost::numeric::ublas::vector container

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

template <class T, class STORAGE>
inline alps::IDump& operator >> (alps::IDump& dump,
                                   boost::numeric::ublas::vector<T,STORAGE>& x)
{
  x.resize(uint32_t(dump));
  if (x.size())
    dump.read_array(x.size(),&(x[0]));
  return dump;
}

/// serialize a boost::numeric::ublas::vector container
template <class T, class STORAGE>
inline alps::ODump& operator << (alps::ODump& dump,
                                   const boost::numeric::ublas::vector<T,STORAGE>& x)
{
  dump << uint32_t(x.size());
  if(x.size())
    dump.write_array(x.size(),&(x[0]));
  return dump;
}          

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_BOOST_UBLAS_H
