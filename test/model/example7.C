/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Axel Grzesik <axel@th.physik.uni-bonn.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/model.h>
#include <fstream>
#include <iostream>
#include <string>

void write_two_site_basis(const std::string& name, const alps::ModelLibrary& lib, 
                          const alps::Parameters& p=alps::Parameters())
{
  alps::SiteBasisDescriptor<short> sitebasis=lib.get_site_basis(name),
                                   sitebasis2=lib.get_site_basis(name);
  sitebasis.set_parameters(p);
  std::cout << "States of basis " << name << "=" << alps::site_basis<short>(sitebasis);
  for(int i=0;i<sitebasis.size();++i)
    sitebasis[i]+=sitebasis2[i];
  std::cout << "States of basis " << name << "=" << alps::site_basis<short>(sitebasis);
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  
  std::ifstream in("../../lib/xml/models.xml");
  alps::ModelLibrary lib(in);
  alps::Parameters p;
  p["NMax"]=2;
  write_two_site_basis("spinful boson",lib,p);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  //alps::comm_exit(true);
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
