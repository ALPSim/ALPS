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
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
