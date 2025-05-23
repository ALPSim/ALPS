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

#include <alps/parapack/clone_info.h>
#include <alps/parapack/clone_info_p.h>
#include <boost/filesystem/operations.hpp>
#include <iostream>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::clone_info info;
  alps::clone_info_xml_handler handler(info);

  alps::XMLParser parser(handler);
  parser.parse(std::cin);

  alps::oxstream ox(std::cout);

  ox << info;

  boost::filesystem::path xdrpath("clone_info.xdr");
  {
    alps::OXDRFileDump dp(xdrpath);
    dp << info;
  }
  info = alps::clone_info();
  {
    alps::IXDRFileDump dp(xdrpath);
    dp >> info;
  }
  ox << info;
  boost::filesystem::remove(xdrpath);

  boost::filesystem::path h5path("clone_info.h5");
  #pragma omp critical (hdf5io)
  {
    alps::hdf5::archive ar(h5path.string(), "a");
    ar["/info"] << info;
  }
  info = alps::clone_info();
  #pragma omp critical (hdf5io)
  {
    alps::hdf5::archive ar(h5path.string());
    ar["/info"] >> info;
  }
  ox << info;
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
