/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                           Matthias Troyer <troyer@comp-phys.org>                *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_HPP
#define ALPS_NGS_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/hdf5/pair.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/pointer.hpp>
#include <alps/hdf5/complex.hpp>

#include <alps/ngs/api.hpp>
#include <alps/ngs/cast.hpp>
#include <alps/ngs/sleep.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/random01.hpp>
#include <alps/ngs/boost_mpi.hpp>
#include <alps/ngs/short_print.hpp>
#include <alps/ngs/thread_exceptions.hpp>
#include <alps/ngs/observablewrappers.hpp> // TODO: remove!

#ifdef ALPS_NGS_USE_NEW_ALEA
	#include <alps/ngs/accumulator.hpp>
#else
	namespace alps {
		namespace accumulator {

			typedef ::alps::ngs::SimpleRealObservable SimpleRealObservable;
			typedef ::alps::ngs::SimpleRealVectorObservable SimpleRealVectorObservable;

			typedef ::alps::ngs::RealObservable RealObservable;
			typedef ::alps::ngs::RealVectorObservable RealVectorObservable;

			typedef ::alps::ngs::SignedRealObservable SignedRealObservable;
			typedef ::alps::ngs::SignedRealVectorObservable SignedRealVectorObservable;

			typedef ::alps::ngs::SignedSimpleRealObservable SignedSimpleRealObservable;
			typedef ::alps::ngs::SignedSimpleRealVectorObservable SignedSimpleRealVectorObservable;
		}
	}
#endif

// #include <alps/mcbase.hpp>
// #include <alps/parseargs.hpp>
// #include <alps/stop_callback.hpp>
// #include <alps/progress_callback.hpp> // TODO: remove this file!

// TODO: remove these deprecated headers:
#include <alps/ngs/mcresult.hpp>
#include <alps/ngs/mcresults.hpp>
#include <alps/ngs/mcoptions.hpp>
#include <alps/ngs/mcobservable.hpp>
#include <alps/ngs/mcobservables.hpp> // TODO: rethink this!

#endif
