/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "fulldiag.h"
#include <alps/utility/copyright.hpp>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

namespace {

void print_usage(std::ostream& out, const char* pname)
{
  out << "Usage:\n"
      << pname << " [--T_MIN ...] [--T_MAX ...] [--DELTA_T ...] [--H_MIN ...] [--H_MAX ... ] [--DELTA_H ... ] [--versus h] [--DENSITIES ...] filenames\n"
      << "or:\n"
      << pname << " --couple mu [--T_MIN ...] [--T_MAX ...] [--DELTA_T ...] [--MU_MIN ...] [--MU_MAX ... ] [--DELTA_MU ...] [--versus mu] [--DENSITIES ...] filenames\n"
      << "\nOptions:\n"
      << "  -h, --help     produce help message\n"
      << "  -l, --license  print license conditions\n";
}

} // namespace

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  int i=1;  
  alps::Parameters parms;

  while (i<argc && argv[i][0]=='-') {
    if (!std::strcmp(argv[i], "--help") || !std::strcmp(argv[i], "-h")) {
      print_usage(std::cout, argv[0]);
      alps::print_copyright(std::cout);
      return 0;
    }
    if (!std::strcmp(argv[i], "--license") || !std::strcmp(argv[i], "-l")) {
      alps::print_license(std::cout);
      return 0;
    }
    if (!std::strcmp(argv[i], "--")) {
      ++i;
      break;
    }
    if (argv[i][1]!='-' || argv[i][2]=='\0') {
      std::cerr << "Unknown option: " << argv[i] << "\n";
      print_usage(std::cerr, argv[0]);
      return 1;
    }
    if (i+1>=argc) {
      std::cerr << "Missing value for option: " << argv[i] << "\n";
      print_usage(std::cerr, argv[0]);
      return 1;
    }
    parms[argv[i]+2]=argv[i+1];
    i+=2;
  }

  alps::print_copyright(std::cout);

  // no filename found
  if(i >= argc) {
    print_usage(std::cerr, argv[0]);
    return 1;
  }

  while (i<argc) {
    boost::filesystem::path p(argv[i]);
    std::string name=argv[i];
    name.erase(name.rfind(".out.xml"),8);
    alps::ProcessList nowhere;
    FullDiagMatrix<double> matrix (nowhere,p);
    matrix.evaluate(parms,name); 
    ++i; 
  }
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif
}
