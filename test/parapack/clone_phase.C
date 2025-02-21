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
