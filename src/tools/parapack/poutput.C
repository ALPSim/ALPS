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

#include <alps/parapack/clone.h>
#include <alps/scheduler/montecarlo.h>

int main(int argc, char** argv) {
  alps::oxstream os(std::cout);
  os << alps::start_tag("CLONES");
  for (int i = 1; i < argc; ++i) {
    alps::hdf5::archive ar(argv[i]);
    alps::Parameters params;
    alps::clone_info info;
    ar["/parameters"] >> params;
    ar["/log/alps"] >> info;
    std::vector<alps::alea::observable_set> obs;
    alps::load_observable(ar, info.clone_id(), obs);
    os << alps::start_tag("CLONE")
       << alps::attribute("dumpfile", argv[i])
       << params;
    for (alps::alea::observable_set const& m : obs) alps::scheduler::write_xml_averages(os, m.results());
    info.write_xml(os);
    os << alps::end_tag("CLONE");
  }
  os << alps::end_tag("CLONES");
}
