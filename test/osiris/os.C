/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/utility/os.hpp>
#include <alps/version.h>
#include <iostream>

int main() {
  std::cerr << "cofig hostname = " << ALPS_CONFIG_HOST << std::endl
            << "config username = " << ALPS_CONFIG_USER << std::endl
            << "exec hostname = " << alps::hostname() << std::endl
            << "exec username = " << alps::username() << std::endl;
  return 0;
}
