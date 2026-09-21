/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/scheduler/convert.h>
#include <alps/parser/xslt_path.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <stdexcept>
#include <iostream>
#include <set>


int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  if (argc<2) {
    std::cerr << "Usage: " << argv[0] << " inputfile [inputfile ...]]\n";
    std::exit(-1);
  }

  std::set<boost::filesystem::path> paths;
  
  for (int i=1;i<argc;++i) {
    std::string inname=argv[i];
    if (inname.size() >= 2 && inname.substr(0, 2) == "./") 
      inname.erase(0, 2);
    alps::convert2xml(inname);
    boost::filesystem::path dir = boost::filesystem::path(inname).remove_filename();
    if (paths.find(dir)==paths.end()) {
      alps::copy_stylesheet(dir);
      paths.insert(dir);
    }
  }

  // make sure ths stylesheet is there

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif

}
