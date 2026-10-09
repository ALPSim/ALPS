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

#include <alps/alea.h>
#include <boost/mpi.hpp>
#include <boost/random.hpp>
#include <iostream>

namespace mpi = boost::mpi;

int main(int argc, char** argv) {
  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.size() >= 2) {
    if (world.rank() == 0) {
      boost::mt19937 engine(2873u);
      boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
        random(engine, boost::uniform_real<>());

      // observables
      alps::ObservableSet obs;
      obs << alps::RealObservable("observable a");
      obs << alps::RealObservable("observable b");

      for(int i=0; i < (1<<12); ++i) {
        obs["observable a"] << random();
        obs["observable b"] << random();
      }
      alps::RealObsevaluator eval_a = obs["observable a"];
      alps::RealObsevaluator eval_b = obs["observable b"];
      alps::RealObsevaluator ratio = eval_a / eval_b;
      obs.addObservable(ratio);
      std::cout << "[output from rank 0]\n" << obs << std::flush;
      world.barrier();

      // send obs to rank 1
      world.send(1, 0, obs);
    } else if (world.rank() == 1) {
      alps::ObservableSet obs;
      world.barrier();

      // receive obs from rank 0
      world.recv(0, 0, obs);
      std::cout << "[output from rank 1]\n" << obs;
    } else {
      // nothing to do
      world.barrier();
    }
  }
}
