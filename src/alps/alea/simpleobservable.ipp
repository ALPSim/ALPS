/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#ifndef ALPS_ALEA_SIMPLEOBSERVABLE_IPP
#define ALPS_ALEA_SIMPLEOBSERVABLE_IPP

#include <alps/alea/simpleobservable.h>

namespace alps {

template <class T,class BINNING> 
void SimpleObservable<T,BINNING>::save(hdf5::archive & ar) const 
{
  AbstractSimpleObservable<T>::save(ar);
  ar[""] << b_;
}

template <class T,class BINNING> 
void SimpleObservable<T,BINNING>::load(hdf5::archive & ar) 
{
  AbstractSimpleObservable<T>::load(ar);
  ar[""] >> b_;
}

template <class T,class BINNING> 
hdf5::archive & operator<<(hdf5::archive & ar,  SimpleObservable<T,BINNING> const& obs) 
{
  return ar["/simulation/results/" + obs.representation()] << obs;
}

template <class T,class BINNING> 
hdf5::archive & operator>>(hdf5::archive & ar,  SimpleObservable<T,BINNING>& obs) 
{
  return ar["/simulation/results/" + obs.representation()] >> obs;
}

} // end namespace alps

#endif // ALPS_ALEA_SIMPLEOBSERVABLE_IPP
