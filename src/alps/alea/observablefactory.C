/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

/* $Id: observableset.C 3765 2010-01-22 14:58:52Z troyer $ */

#include <alps/config.h>

#include "observablefactory.h"
#include "signedobservable.h"
#include "nobinning.h"
#include "detailedbinning.h"
#include "simpleobseval.h"
#include "histogrameval.h"

#include "simpleobseval.ipp"
#include "abstractsimpleobservable.ipp"
#include "simpleobservable.ipp"

namespace alps {

// clang > 3.3 does not find the observable typedef, so instanciat it here ...  
template class SimpleObservable<int32_t,DetailedBinning<int32_t> >; // IntObservable;
template class SimpleObservable<double,DetailedBinning<double> >; // RealObservable;
// template class SimpleObservable<float,DetailedBinning<float> >; // FloatObservable;
// template class SimpleObservable<std::complex<double>,DetailedBinning<std::complex<double> > >; // ComplexObservable;
template class SimpleObservable<double,FixedBinning<double> >; // RealTimeSeriesObservable;
template class SimpleObservable<int32_t,FixedBinning<int32_t> >; // IntTimeSeriesObservable;
template class SimpleObservable< std::valarray<int32_t> , 
                         DetailedBinning<std::valarray<int32_t> > >; // IntVectorObservable;
template class SimpleObservable< std::valarray<double> , 
                         DetailedBinning<std::valarray<double> > >; // RealVectorObservable;
// template class SimpleObservable< std::valarray<float> , 
                         // DetailedBinning<std::valarray<float> > >; // FloatVectorObservable;
//template class SimpleObservable< std::valarray<std::complex<double> > , 
//                         DetailedBinning<std::valarray<std::complex<double> > > >; // ComplexVectorObservable;
template class SimpleObservable< std::valarray<int32_t> , 
                         FixedBinning<std::valarray<int32_t> > >; // IntVectorTimeSeriesObservable;
template class SimpleObservable< std::valarray<double> , 
                         FixedBinning<std::valarray<double> > >; // RealVectorTimeSeriesObservable;
//template class SimpleObservable< std::valarray<std::complex<double> > , 
//                         FixedBinning<std::valarray<std::complex<double> > > >; // ComplexVectorTimeSeriesObservable;

template class SimpleObservable<int32_t,NoBinning<int32_t> >; // SimpleIntObservable;
template class SimpleObservable<double,NoBinning<double> >; // SimpleRealObservable;
// template class SimpleObservable<float,NoBinning<float> >; // SimpleFloatObservable;
// template class SimpleObservable<std::complex<double>,NoBinning<std::complex<double> > >; // SimpleComplexObservable;
template class SimpleObservable< std::valarray<int32_t> , NoBinning<std::valarray<int32_t> > >; // SimpleIntVectorObservable;
template class SimpleObservable< std::valarray<double> , NoBinning<std::valarray<double> > >; // SimpleRealVectorObservable;
// template class SimpleObservable< std::valarray<float> , NoBinning<std::valarray<float> > >; // SimpleFloatVectorObservable;
// template class SimpleObservable< std::valarray<std::complex<double> > ,
//                          NoBinning<std::valarray<std::complex<double> > > >; // SimpleComplexVectorObservable;

// instanciate the base classes of the observables
template class AbstractSimpleObservable<double>;
template class AbstractSimpleObservable<std::valarray<double> >;
// instanciate evaluators
template class SimpleObservableEvaluator<double>;
template class SimpleObservableEvaluator<std::valarray<double> >;


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
