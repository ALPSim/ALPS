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

#include <alps/parapack/clone_info_p.h>
#include <alps/parser/xmlparser.h>
#include <boost/filesystem/operations.hpp>
#include <iostream>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::clone_phase phase;
  alps::clone_phase_xml_handler handler(phase);

  alps::XMLParser parser(handler);
  parser.parse(std::cin);

  alps::oxstream ox(std::cout);

  ox << phase;

  boost::filesystem::path xdrpath("clone_phase.xdr");
  {
    alps::OXDRFileDump dp(xdrpath);
    dp << phase;
  }
  phase = alps::clone_phase();
  {
    alps::IXDRFileDump dp(xdrpath);
    dp >> phase;
  }
  ox << phase;
  boost::filesystem::remove(xdrpath);

  boost::filesystem::path h5path("clone_phase.h5");
  #pragma omp critical (hdf5io)
  {
    alps::hdf5::archive ar(h5path.string(), "a");
    ar["/phase"] << phase;
  }
  phase = alps::clone_phase();
  #pragma omp critical (hdf5io)
  {
    alps::hdf5::archive ar(h5path.string());
    ar["/phase"] >> phase;
  }
  ox << phase;
  boost::filesystem::remove(h5path);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
