/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2009 by Synge Todo <wistaria@comp-phys.org>
*
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <boost/version.hpp>
#if BOOST_VERSION >= 103600
# include <boost/detail/workaround.hpp>
# if !defined(BOOST_SPIRIT_USE_OLD_NAMESPACE)
#  define BOOST_SPIRIT_USE_OLD_NAMESPACE
# endif
# include <boost/spirit/include/classic_actor.hpp>
# include <boost/spirit/include/classic_core.hpp>
# include <boost/spirit/include/classic_confix.hpp>
#else
# include <boost/spirit/actor.hpp>
# include <boost/spirit/core.hpp>
# include <boost/spirit/utility/confix.hpp>
#endif
