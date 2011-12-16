/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
