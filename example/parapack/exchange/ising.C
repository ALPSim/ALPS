/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#include "../single/ising.h"
#include <alps/parapack/exchange.h>
#ifdef ALPS_HAVE_MPI
#include "../multiple/ising.h"
#include <alps/parapack/exchange_multi.h>
#endif

PARAPACK_SET_VERSION("ALPS/parapack example program: exchange Monte Carlo");
PARAPACK_REGISTER_ALGORITHM(single_ising_worker, "ising");
PARAPACK_REGISTER_EVALUATOR(ising_evaluator, "ising");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<single_ising_worker>,
                         "ising; exchange");
#ifdef ALPS_HAVE_MPI
PARAPACK_REGISTER_PARALLEL_WORKER(alps::parapack::parallel_exchange_worker<single_ising_worker>,
                                  "ising; exchange");
PARAPACK_REGISTER_PARALLEL_WORKER(alps::parapack::multiple_parallel_exchange_worker<parallel_ising_worker>,
                                  "multiple parallel ising; exchange");
#endif
PARAPACK_REGISTER_EVALUATOR(ising_evaluator, "ising; exchange");
