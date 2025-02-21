/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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
