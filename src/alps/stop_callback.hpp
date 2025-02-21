/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>,                  *
 *                              Synge Todo <wistaria@comp-phys.org>                *
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

#ifndef ALPS_NGS_CALLBACK_HPP
#define ALPS_NGS_CALLBACK_HPP

#include <alps/ngs/config.hpp>
#include <alps/ngs/signal.hpp>

#include <boost/chrono.hpp>
#ifdef ALPS_HAVE_MPI
# include <boost/mpi/communicator.hpp>
#endif

namespace alps {

	class ALPS_DECL stop_callback {
		public:
		    stop_callback(std::size_t timelimit);
#ifdef ALPS_HAVE_MPI
			stop_callback(boost::mpi::communicator const & cm, std::size_t timelimit);
#endif
		    bool operator()();
		private:
		    boost::chrono::duration<std::size_t> limit;
		    alps::ngs::signal signals;
		    boost::chrono::high_resolution_clock::time_point start;
#ifdef ALPS_HAVE_MPI
	        boost::optional<boost::mpi::communicator> comm;
#endif
	};

#ifdef ALPS_HAVE_MPI
		// TODO: remove this!
        class ALPS_DECL stop_callback_mpi {
        public:
          stop_callback_mpi(boost::mpi::communicator const& cm, std::size_t timelimit);
          bool operator()();
        private:
          boost::mpi::communicator comm;
          boost::chrono::duration<std::size_t> limit;
          alps::ngs::signal signals;
          boost::chrono::high_resolution_clock::time_point start;
	};
#endif
}

#endif
