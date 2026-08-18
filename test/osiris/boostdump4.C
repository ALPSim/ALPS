/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
    alps::odump_archive ar(od);
    ar << o1 << o2 << o3 << o4 << o5 << o6 << o7 << o8 << o9 << o10 << o11;
  }
  
  alps::IXDRFileDump id(boost::filesystem::path("xdrdump.dump"));
  std::cout << id.get<bool>() << ' ';
  std::cout << static_cast<int32_t>(id.get<int8_t>()) << ' ';
  std::cout << static_cast<int32_t>(id.get<uint8_t>()) << ' ';
  std::cout << id.get<int16_t>() << ' ';
  std::cout << id.get<uint16_t>() << ' ';
  std::cout << static_cast<int32_t>(id) << ' ';
  std::cout << static_cast<uint32_t>(id) << ' ';
  int64_t i8 = id;
  uint64_t i9(id);
  std::cout << i8 << ' '  << i9  << ' ';
  double i10 = static_cast<double>(id);
  std::cout << i10 << ' ';
  std::string str;
  id >> str;
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
