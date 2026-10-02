/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/clone_info.h>
#include <alps/osiris/comm.h>
#include <chrono>
#include <thread>

int main(int argc, char **argv) {
  alps::comm_init(argc, argv);
  alps::Parameters params;
  params["SEED"] = 29832;
  alps::clone_info info(0, params, "info_test");
  info.start("test 1");
  std::this_thread::sleep_for(std::chrono::seconds(1));
  info.stop();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  info.start("test 2");
  std::this_thread::sleep_for(std::chrono::seconds(1));
  info.stop();
  info.set_progress(0.593483);
  if (alps::is_master()) {
    alps::oxstream oxs;
    oxs << info;
  }
  info.set_progress(1);
  if (alps::is_master()) {
    alps::oxstream oxs;
    oxs << info;
  }
  alps::comm_exit();
  return 0;
}
