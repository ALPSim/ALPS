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

  int np = world.size();
  int pid = world.rank();

  if (np < 2) std::exit(-1);

  boost::mt19937 engine;

  std::vector<int> nelms, offsets;
  int n = 0;
  for (int p = 0; p < np; ++p) {
    nelms.push_back((engine() & 3) + 1);
    offsets.push_back(n);
    n += nelms.back();
  }

  std::vector<uint32_t> source;
  if (world.rank() == 0) {
    for (int i = 0; i < n; ++i) source.push_back(engine());
    std::cout << "source: " << alps::write_vector(source) << std::endl;
  }

  std::vector<uint32_t> distributed;
  alps::distribute_vector(world, nelms, offsets, source, distributed);
  for (int p = 0; p < np; ++p) {
    if (p == pid)
      std::cout << "process " << p << ": " << alps::write_vector(distributed) << std::endl;
    std::cout << std::flush;
    world.barrier();
  }

  std::vector<uint32_t> collected(world.rank() == 0 ? n : 0);
  alps::collect_vector(world, nelms, offsets, distributed, collected);
  if (world.rank() == 0)
    std::cout << "collected: " << alps::write_vector(collected) << std::endl;

  if (world.rank() == 0 && source != collected) {
    std::cout << "test failed\n";
    std::exit(-1);
  }
}
