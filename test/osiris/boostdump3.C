/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#include <alps/osiris.h>
#include <iostream>
#include <cstdlib>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  bool o1 = false;
  int8_t o2 = 63;
  uint8_t o3 = 201;
  int16_t o4 = -699;
  uint16_t o5 = 43299;
  int32_t o6 = 847229;
  uint32_t o7 = 4294967295u;
  int64_t o8 = -1152921504606846976ll;
  uint64_t o9 = 18446744073709551614ull;
  double o10 = 3.14159265358979323846;
  std::string o11 = "test string";

  {
    alps::OXDRFileDump od(boost::filesystem::path("xdrdump.dump"));
    od << o1 << o2 << o3 << o4 << o5 << o6 << o7 << o8 << o9 << o10 << o11;
  }
  
  alps::IXDRFileDump id(boost::filesystem::path("xdrdump.dump"));
  alps::idump_archive ar(id);
  std::cout << alps::get<bool>(ar) << ' ';
  std::cout << static_cast<int32_t>(alps::get<int8_t>(ar)) << ' ';
  std::cout << static_cast<int32_t>(alps::get<uint8_t>(ar)) << ' ';
  std::cout << alps::get<int16_t>(ar) << ' ';
  std::cout << alps::get<uint16_t>(ar) << ' ';
  std::cout << alps::get<int32_t>(ar) << ' ';
  std::cout << alps::get<uint32_t>(ar) << ' ';
  int64_t i8 = alps::get<int64_t>(ar);
  uint64_t i9(alps::get<uint64_t>(ar));
  std::cout << i8 << ' '  << i9  << ' ';
  double i10 = alps::get<double>(ar);
  std::cout << i10 << ' ';
  std::string str;
  ar >> str;
  std::cout << str << std::endl;
  
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
