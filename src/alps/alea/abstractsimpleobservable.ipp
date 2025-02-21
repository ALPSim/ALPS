/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2012 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@itp.phys.ethz.ch>,
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

#ifndef ALPS_ALEA_ABSTRACTSIMPLEOBSERVABLE_IPP
#define ALPS_ALEA_ABSTRACTSIMPLEOBSERVABLE_IPP

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/valarray.hpp>

#include <alps/alea/abstractsimpleobservable.h>

namespace alps {

template <class T>
void AbstractSimpleObservable<T>::save(hdf5::archive & ar) const 
{
  Observable::save(ar);
  if (label_.size())
      ar
          << make_pvp("labels", label_)
      ;
  ar
      << make_pvp("count", count())
      ;
  if (count() > 0) {
      ar
          << make_pvp("mean/value", mean())
      ;
  }
  if (count() > 1) {
      ar
          << make_pvp("mean/error", error())
          << make_pvp("mean/error_convergence", converged_errors())
      ;
      if(has_variance())
          ar
              << make_pvp("variance/value", variance())
          ;
      if(has_tau())
          ar
              << make_pvp("tau/value", tau())
          ;
  }
}

template <class T>
void AbstractSimpleObservable<T>::load(hdf5::archive & ar) 
{
    Observable::load(ar);
    if (ar.is_data("labels"))
        ar >> make_pvp("labels", label_);
}

}

#endif // ALPS_ALEA_ABSTRACTSIMPLEOBSERVABLE_IPP
