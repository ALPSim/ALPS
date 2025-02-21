/* 
 * Copyright Matthias Troyer 2006
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
 */


#ifndef ALPS_RANDOM_PARALLEL_MPI_HPP
#define ALPS_RANDOM_PARALLEL_MPI_HPP

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <alps/random/parallel/seed.hpp>

/// @file This header provides convenience seeding functions for applications based on the candidate 
/// Boost.MPI library. The @c num and @c total parameters are obtained from a 
/// @c boost::parallel::mpi::communicator MPI communicator.

namespace alps { namespace random { namespace parallel {
  
  /// This seed function calls the corresponding @c seed functions in the header @c alps/parallel/seed.hpp
  /// using the communicator rank @c c.rank() and size @c c.size()
  /// as the @c num  and @c total arguments.

  template <class PRNG>
  void seed(
          PRNG& prng
        , boost::mpi::communicator const& c
      )
  {
    seed(prng, c.rank(), c.size());
  }

  /// This seed function calls the corresponding @c seed functions in the header @c alps/parallel/seed.hpp
  /// using the communicator rank @c c.rank() and size @c c.size()
  /// as the @c num  and @c total arguments.
  template <class PRNG, class SeedType>
  void seed(
          PRNG& prng
        , boost::mpi::communicator const& c
        , SeedType const& s
      )
  {
    seed(prng, c.rank(), c.size(), s);
  }

  /// This seed function calls the corresponding @c seed functions in the header @c alps/parallel/seed.hpp
  /// using the communicator rank @c c.rank() and size @c c.size()
  /// as the @c num  and @c total arguments.
  template <class PRNG, class Iterator>
  void seed(
          PRNG& prng
        , boost::mpi::communicator const& c
        , Iterator& first
        , Iterator const& last
      )
  {
    seed(prng, c.rank(), c.size(), first, last);
  }

/// broadcasts the global seed @c s from the root rank given by the @c root argument and then
/// calls the seed function with the same arguments.

  template <class PRNG, class SeedType>
  void broadcast_seed(
          PRNG& prng
        , boost::mpi::communicator const& c
        , int root
        , SeedType s
      )
  {
    boost::mpi::broadcast(c,s,root);
    seed(prng, c.rank(), c.size(), s);
  }

} } } // namespace alps::random::parallel

#endif // ALPS_RANDOM_PARALLEL_MPI_HPP
