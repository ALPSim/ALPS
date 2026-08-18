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
  mpi::environment env(argc, argv);
  mpi::communicator world;

  alps::process_helper_mpi process(world, 4);

  if (world.rank() == 0) {
    alps::process_group g1 = process.allocate();
    for (int i = 0; i < g1.process_list.size(); ++i)
      std::cout << g1.process_list[i] << ' ';
    std::cout << std::endl;
    std::cout << process.num_groups() << ' ' << process.num_free() << std::endl;

    alps::process_group g2 = process.allocate();
    for (int i = 0; i < g2.process_list.size(); ++i)
      std::cout << g2.process_list[i] << ' ';
    std::cout << std::endl;
    std::cout << process.num_groups() << ' ' << process.num_free() << std::endl;

    process.release(g1);
    std::cout << process.num_groups() << ' ' << process.num_free() << std::endl;

    g1 = process.allocate();
    for (int i = 0; i < g1.process_list.size(); ++i)
      std::cout << g1.process_list[i] << ' ';
    std::cout << std::endl;
    std::cout << process.num_groups() << ' ' << process.num_free() << std::endl;

    process.release(g2);
    std::cout << process.num_groups() << ' ' << process.num_free() << std::endl;

    process.release(g1);
    std::cout << process.num_groups() << ' ' << process.num_free() << std::endl;
  }

  process.halt();
  while (true) {
    if (process.check_halted()) break;
  }
}
