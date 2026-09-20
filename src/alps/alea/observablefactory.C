/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: observableset.C 3765 2010-01-22 14:58:52Z troyer $ */

#include <alps/config.h>

#include "observablefactory.h"
#include "signedobservable.h"
#include "nobinning.h"
#include "detailedbinning.h"
#include "simpleobseval.h"
#include "histogrameval.h"


namespace alps {

ObservableFactory::ObservableFactory()
{
  register_observable<IntObsevaluator>();
  register_observable<RealObsevaluator>();
  register_observable<IntObservable>();
  register_observable<RealObservable>();
  register_observable<AbstractSignedObservable<RealObsevaluator> >();
  register_observable<SignedObservable<RealObservable> >();
  register_observable<SignedObservable<SimpleRealObservable> >();
  register_observable<SignedObservable<RealTimeSeriesObservable> >();
  register_observable<IntTimeSeriesObservable>();
  register_observable<RealTimeSeriesObservable>();
  register_observable<SimpleRealObservable>();
  register_observable<SimpleIntObservable>();
  register_observable<RealVectorObsevaluator>();
  register_observable<RealVectorObservable>();
  register_observable<SignedObservable<RealVectorObservable> >();
  register_observable<IntVectorObservable>();
  register_observable<IntVectorObsevaluator>();
  register_observable<IntVectorTimeSeriesObservable>();
  register_observable<RealVectorTimeSeriesObservable>();
  register_observable<SimpleIntVectorObservable>();
  register_observable<SimpleRealVectorObservable>();
  register_observable<AbstractSignedObservable<RealVectorObsevaluator> >();
  register_observable<SignedObservable<SimpleRealVectorObservable> >();
  register_observable<SignedObservable<RealVectorTimeSeriesObservable> >();
  register_observable<IntHistogramObservable>();
  register_observable<RealHistogramObservable>();
  register_observable<IntHistogramObsevaluator>();
  register_observable<RealHistogramObsevaluator>();
/*
  register_observable<HistogramObservable<int32_t> >();
  register_observable<HistogramObservable<int32_t,double> >();
*/
}

} // namespace alps
