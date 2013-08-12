/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
