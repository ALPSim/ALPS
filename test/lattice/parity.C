/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#include <alps/lattice.h>
#include <alps/parameter.h>
#include <fstream>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main() {
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    typedef alps::coordinate_graph_type graph_t;
    typedef alps::parity_t parity_t;

    // read parameters
    alps::ParameterList parms;
    std::cin >> parms;
    
    for (alps::ParameterList::const_iterator p = parms.begin();
         p != parms.end(); ++p) {
      // create the lattice

      alps::graph_helper<> lattice(*p);
      const graph_t& graph = lattice.graph();
      
      std::cout << graph;
      
      for (graph_t::vertex_iterator vi = boost::vertices(graph).first;
           vi != boost::vertices(graph).second; ++vi) {
        std::cout << "vertex " << *vi << "'s parity is ";
        if (boost::get(parity_t(), graph, *vi)
            == alps::parity_traits<parity_t, graph_t>::white) {
          std::cout << "white\n";
        } else if (boost::get(alps::parity_t(), graph, *vi)
                   == alps::parity_traits<parity_t, graph_t>::black) {
          std::cout << "black\n";
        } else {
          std::cout << "undefined\n";
        }        
      }
    }
         
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    exit(-1);
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    exit(-2);
  }
#endif
  return 0;
}
