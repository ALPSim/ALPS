/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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
#include <alps/osiris/archivedump.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <iostream>
#include <cstdlib>
#include <fstream>

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
    std::ofstream of("textdump.dump");
	if (!of)
		std::cerr << "Issue\n";
	boost::archive::text_oarchive oa(of);//of,boost::archive::no_codecvt);
    alps::archive_odump<boost::archive::text_oarchive> od(oa);
    od << o1 << o2 << o3 << o4 << o5 << o6 << o7 << o8 << o9 << o10 << o11;
  }

  std::ifstream inf("textdump.dump");
  boost::archive::text_iarchive ia(inf);
  alps::archive_idump<boost::archive::text_iarchive> id(ia);

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
