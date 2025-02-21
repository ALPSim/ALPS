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

/* $Id: example4.C 730 2004-04-02 18:20:12Z troyer $ */

#include <alps/lattice.h>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main(int argc, char ** argv)
{

#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    alps::Parameters parameters;
    switch(argc) {
       case 1:
         std::cin >> parameters;
         break;
       case 2:
         {
           std::string fn(argv[1]);
           std::ifstream parmfile(fn.c_str());

           if(boost::algorithm::iends_with(fn, ".xml"))
             parameters.extract_from_xml(parmfile);
           else
             parmfile >> parameters;
         }
         break;
       default:
         std::cerr << "Usage: " << argv[0] << " [parameterfile]\n";
    } 
   
    // create a graph factory with default graph type
    alps::graph_helper<> lattice(parameters);
    // write the graph created from the input in XML
    std::cout << lattice.graph();

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
}
