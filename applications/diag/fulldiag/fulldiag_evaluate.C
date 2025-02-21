/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
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

#include "fulldiag.h"
#include <fstream>

// Print usage and exit

void error_exit(char *pname)
 {
    std::cerr << "Usage:\n" << pname << " [--T_MIN ...] [--T_MAX ...] [--DELTA_T ...] [--H_MIN ...] [--H_MAX ... ] [--DELTA_H ... ] [--versus h] [--DENSITIES ...] filenames\n";
    std::cerr << "or:\n" << pname << " --couple mu [--T_MIN ...] [--T_MAX ...] [--DELTA_T ...] [--MU_MIN ...] [--MU_MAX ... ] [--DELTA_MU ...] [--versus mu] [--DENSITIES ...] filenames\n";
    exit(1);
 }

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  int i=1;  
  alps::Parameters parms;

  while (i<argc-1 && argv[i][0]=='-') {
    parms[argv[i]+2]=argv[i+1];
    i+=2;
  }

  // no filename found
  if(i >= argc)
    error_exit(argv[0]);

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
