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

#include <alps/parapack/clone_info.h>
#include <alps/osiris/comm.h>

namespace mpi = boost::mpi;

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  mpi::communicator world;
  alps::Parameters params;
  params["SEED"] = 29832;
  alps::clone_info_mpi info(world, 0, params, "info_test");
  info.start("test 1");
  sleep(1);
  info.stop();
  sleep(1);
  info.start("test 2");
  sleep(1);
  info.stop();
  info.set_progress(0.593483);
  if (world.rank() == 0) {
    alps::oxstream oxs;
    oxs << info;
  }
  info.set_progress(1);
  if (world.rank() == 0) {
    alps::oxstream oxs;
    oxs << info;
  }
}
