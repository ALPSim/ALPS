/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/process.h>
#include <iostream>

namespace mpi = boost::mpi;

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  mpi::communicator world;
  if (world.size() >= 4) {
    alps::process_helper_mpi process(world, 4);
    mpi::communicator cg = process.comm_ctrl();
    mpi::communicator cl = process.comm_work();
    mpi::communicator hd = process.comm_head();
    for (int p = 0; p < world.size(); ++p) {
      if (world.rank() == p) {
        std::cout << "rank: " << world.rank();
        if (cg)
          std::cout << ", global rank = " << cg.rank();
        if (cl)
          std::cout << ", local rank = " << cl.rank();
        if (hd)
          std::cout << ", head rank = " << hd.rank();
        std::cout << std::endl;
      }
      std::cout << std::flush;
      world.barrier();
    }
    process.halt();
    while (!process.check_halted()) {}
  }
}
