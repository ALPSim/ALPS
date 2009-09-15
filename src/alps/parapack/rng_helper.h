/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef PARAPACK_RNG_HELPER_H
#define PARAPACK_RNG_HELPER_H

#include "types.h"
#include <alps/config.h>
#include <alps/osiris.h>
#include <alps/parameter.h>
#include <alps/random/buffered_rng.h>

namespace alps {

class rng_helper {
public:
  rng_helper(const Parameters& p);
  void init(const Parameters& p);
  void load(IDump& dp);
  void save(ODump& dp) const;
  typedef buffered_rng_base engine_type;
  typedef boost::variate_generator<engine_type&, boost::uniform_real<> > generator_type;
  uint32_t seed;
  uint32_t disorder_seed; // shared by all workers in each clone
  engine_type& engine() { return *(engines_[0]); }
  engine_type& engine(int p) { return *(engines_[p]); }
  generator_type& generator_01() { return (*generators_[0]); }
  generator_type& generator_01(int p) { return (*generators_[p]); }
  double random_01() { return generators_[0]->operator()(); }
  double random_01(int p) { return generators_[p]->operator()(); }
  double uniform_01() { return random_01(); }
  double uniform_01(int p) { return random_01(p); }
  // int random_int(int a, int b) { return a + int((b-a+1) * rngs[0]()); }
  // int random_int(int n) { return int(n * uniform_01()); }
  // double random() { return uniform_01(); } // obsolete
private:
  mutable std::vector<boost::shared_ptr<engine_type> > engines_;
  std::vector<boost::shared_ptr<generator_type> > generators_;
};

} // end namespace alps

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

inline alps::IDump& operator>>(alps::IDump& dp, alps::rng_helper& rng) {
  rng.load(dp);
  return dp;
}

inline alps::ODump& operator<<(alps::ODump& dp, alps::rng_helper const& rng) {
  rng.save(dp);
  return dp;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // PARAPACK_RNG_HELPER_H
