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

/* $Id: observableset.h 3696 2010-01-13 15:25:30Z gamperl $ */

#ifndef ALPS_ALEA_OBSERVABLEFACTORY_H
#define ALPS_ALEA_OBSERVABLEFACTORY_H

#include <alps/config.h>
#include <alps/factory.h>
#include <alps/alea/observable.h>

namespace alps {

/** A class to collect the various measurements performed in a simulation
    It is implemented as a map, with std::string as key type */

class ObservableFactory : public factory<uint32_t,Observable>
{
public:
  ObservableFactory();
  template <class T>
  void register_observable() { register_type<T>(T::version); }
};

} // end namespace alps

#endif // ALPS_ALEA_OBSERVABLEFACTORY_H
