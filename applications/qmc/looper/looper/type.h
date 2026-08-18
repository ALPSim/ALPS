/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef LOOPER_TYPE_H
#define LOOPER_TYPE_H

#include <boost/mpl/bool.hpp>

namespace looper {

//
// QMC types
//

struct classical {};
struct path_integral {};
struct sse {};

//
// meta functions
//

template<typename QMC>
struct is_path_integral
{ typedef boost::mpl::false_ type; };

template<>
struct is_path_integral<path_integral>
{ typedef boost::mpl::true_ type; };

template<typename QMC>
struct is_sse
{ typedef boost::mpl::false_ type; };

template<>
struct is_sse<sse>
{ typedef boost::mpl::true_ type; };

} // end namepspace looper

#endif // LOOPER_TYPE_H
