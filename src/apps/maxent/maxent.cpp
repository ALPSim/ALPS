/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/solvers.hpp>
#include <alps/ngs.hpp>

int main(int argc, char** argv) {
  alps::mcoptions options(argc, argv);
  alps::params parms(alps::hdf5::archive(options.input_file));
  std::string output_file = boost::lexical_cast<std::string>(parms["BASENAME"] | options.output_file) + ".out.h5";
  alps::solvers::maxent(parms, output_file);
  return 0;
}
