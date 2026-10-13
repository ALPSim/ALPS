/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/util.h>
#include <iostream>

int main() {
  std::cout << alps::id2string(0) << std::endl
            << alps::id2string(9) << std::endl
            << alps::id2string(10) << std::endl
            << alps::id2string(15) << std::endl
            << alps::id2string(100) << std::endl
            << alps::id2string(200) << std::endl
            << alps::id2string(1000) << std::endl
            << alps::id2string(1001) << std::endl
            << alps::id2string(100000) << std::endl;
}
