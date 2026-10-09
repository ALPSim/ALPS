/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/random/rngfactory.h>
#include <alps/random/mersenne_twister.hpp>
#include <boost/random.hpp>

alps::RNGFactory::RNGFactory()
{
  register_type<buffered_rng<
    boost::random::lagged_fibonacci<uint32_t, 48, 607, 273>  > >
    ("lagged_fibonacci607");
  register_type<buffered_rng<boost::mt19937> >("mt19937");
}

alps::RNGFactory alps::rng_factory;
