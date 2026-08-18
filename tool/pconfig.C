/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/utility/copyright.hpp>
#include <alps/version.h>
#include <alps/utility/os.hpp>
#include <boost/version.hpp>
#include <iostream>

int main() {
  std::cout << "ALPS version:          " << alps::version() << std::endl
            << "Boost version:         " << BOOST_LIB_VERSION << std::endl
            << "installed at:          " << ALPS_PREFIX << std::endl
            << "configured on:         " << alps::config_host() << std::endl
            << "configured by:         " << alps::config_user() << std::endl
            << "compiled on:           " << alps::compile_date() << std::endl
            << "current hostname:      " << alps::hostname() << std::endl
            << "current user:          " << alps::username() << std::endl;
}
