/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2009 by Synge Todo <wistaria@comp-phys.org>
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

// site coloring test program

#include <alps/lattice.h>
#include <alps/parameter.h>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main() {
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

    typedef alps::graph_helper<> lattice_type;
    typedef lattice_type::graph_type graph_type;
    typedef boost::property_map<graph_type, alps::site_index_t>::const_type vertex_index_map;
    alps::ParameterList params(std::cin);
    BOOST_FOREACH(alps::Parameters const& p, params) {
      lattice_type lattice(p);
      std::vector<std::size_t> color(lattice.num_sites());
      int nc = boost::sequential_vertex_coloring(lattice.graph(),
        boost::iterator_property_map<std::size_t*, vertex_index_map>(&color.front(),
          get(boost::vertex_index, lattice.graph())));
      std::cout << "LATTICE = " << p["LATTICE"] << std::endl;
      std::cout << "  number of colors = " << nc << std::endl;
      std::cout << "  site colors =";
      for (unsigned int s = 0; s < lattice.num_sites(); ++s) std::cout << ' ' << color[s];
      std::cout << std::endl;
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
