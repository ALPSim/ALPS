/* 
 * Copyright Matthias Troyer 2006
 * SPDX-License-Identifier: MIT
*
 */

#ifndef ALPS_RANDOM_PARALLEL_DETAIL_GET_PRIME_HPP
#define ALPS_RANDOM_PARALLEL_DETAIL_GET_PRIME_HPP

#include <boost/cstdint.hpp>
#include <alps/export.h>

namespace alps { namespace random { namespace detail {

// get a prime number to be used as additive constant in a 64-bit LCG generator
    ALPS_DECL boost::uint64_t get_prime_64(unsigned int);

} } } // namespace alps::random::detail

#endif // ALPS_RANDOM_PARALLEL_DETAIL_GET_PRIME_HPP
