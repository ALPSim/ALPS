/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>
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

#include <alps/model.h>
#include <fstream>
#include <iostream>

void write_set(const std::string& name, const alps::ModelLibrary& lib, 
               const alps::Parameters& p=alps::Parameters())
{
  alps::SiteBasisDescriptor<short> sitebasis=lib.get_site_basis(name);
  sitebasis.set_parameters(p);
  std::cout << "States of basis " << name << "=" << alps::site_basis<short>(sitebasis);
}

int main()
{

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    // create the library from an XML file
    std::ifstream in("../../lib/xml/models.xml");
    alps::ModelLibrary lib(in);

    // write all basis states
    write_set("fermion",lib);
    write_set("hardcore boson",lib);
    alps::Parameters p;
    p["local_spin"]=1;
    write_set("spin",lib,p);
    p["local_spin"]=0.5;
    write_set("spin",lib,p);
    p["local_spin"]=1.5;
    write_set("spin",lib,p);
    p["local_spin"]=2;
    write_set("spin",lib,p);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  exit(-1);
}
catch (...)
{
  std::cerr << "Caught unknown exception\n";
  exit(-2);
}
#endif
  return 0;
}
