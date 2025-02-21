/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2010 by Synge Todo <wistaria@comp-phys.org>
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
