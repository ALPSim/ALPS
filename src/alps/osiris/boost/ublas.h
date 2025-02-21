/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2011 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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
