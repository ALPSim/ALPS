/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/process.h>
#include <iostream>

namespace mpi = boost::mpi;

int main(int argc, char** argv) {
  mpi::environment mpi(argc, argv);
  mpi::communicator world;
  alps::process_helper_mpi process(world, 1);
  if (world.rank() == 0) sleep(1);
  process.halt();
  while (true) {
    if (process.check_halted()) {
      std::cerr << "process " << world.rank() << " is halted\n";
      break;
    } else {
      std::cerr << "process " << world.rank() << " is not halted yet\n";
      sleep(1);
    }
  }
}
