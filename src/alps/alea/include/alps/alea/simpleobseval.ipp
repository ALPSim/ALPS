/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_ALEA_SIMPLEOBSEVAL_IPP
#define ALPS_ALEA_SIMPLEOBSEVAL_IPP

#include <alps/alea/simpleobseval.h>
#include <alps/hdf5/valarray.hpp>

namespace alps {

template <typename T> inline void SimpleObservableEvaluator<T>::save(hdf5::archive & ar) const {
    ar[""] << all_;
}
template <typename T> inline void SimpleObservableEvaluator<T>::load(hdf5::archive & ar) {
    ar[""] >> all_;
}

} // end namespace alps

#endif // ALPS_ALEA_SIMPLEOBSEVAL_IPP
