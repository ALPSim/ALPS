/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/lattice.h>
#include <alps/parser/xmlstream.h>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main(int argc, char** argv) {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  // read parameters
  alps::Parameters parameters;
  switch (argc) {
  case 1 :
    std::cin >> parameters;
    break;
  case 2 :
    parameters["LATTICE"] = std::string(argv[1]);
    break;
  default :
    for (int i = 2; i < argc; ++i) {
      std::istringstream iss(argv[i]);
      parameters.parse(iss);
    }
    parameters["LATTICE"] = std::string(argv[1]);
    break;
  }

  // create a graph factory with default graph type
  alps::graph_helper<> lattice(parameters);

  bool unroll = parameters.value_or_default("UNROLL_BOUNDARY", false);

  typedef std::map<std::pair<int /* original id */, int /* crossing direction */>,
    int /* new id */> site_map_t;
  typedef std::vector<boost::tuple<int /* original id */, std::vector<double> /* coordinate */> >
    site_table_t;
  typedef std::vector<boost::tuple<int /* original id */, int /* source */, int /* target */> >
    bond_table_t;

  site_map_t site_map;
  site_table_t outside_sites;
  bond_table_t crossing_bonds;

  if (unroll) {
    int id = lattice.num_sites();
    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()) {
      alps::boundary_crossing crossing = get(alps::boundary_crossing_t(), lattice.graph(), b);
      if (crossing) {
        int s = lattice.source(b);
        int t = lattice.target(b);

        // crossing direction
        int bc_s = 0;
        int bc_t = 0;
        for (std::size_t d = 0; d < lattice.dimension(); ++d) {
          bc_s += (crossing.crosses(d) + 1) << (2*d);
          bc_t += (1- crossing.crosses(d)) << (2*d);
        }

        // coordinate
        alps::graph_helper<>::vector_type cord_s = lattice.coordinate(lattice.target(b));
        alps::graph_helper<>::vector_type cord_t = lattice.coordinate(lattice.source(b));
        alps::graph_helper<>::vector_type vec = lattice.bond_vector(b);
        for (std::size_t d = 0; d < lattice.dimension(); ++d) {
          cord_s[d] -= vec[d];
          cord_t[d] += vec[d];
        }

        // add sites
        int sn;
        if (site_map.find(std::make_pair(s, bc_s)) == site_map.end()) {
          sn = id++;
          site_map[std::make_pair(s, bc_s)] = sn;
          outside_sites.push_back(boost::make_tuple(s, cord_s));
        } else {
          sn = site_map.find(std::make_pair(s, bc_s))->second;
        }

        int tn;
        if (site_map.find(std::make_pair(t, bc_t)) == site_map.end()) {
          tn = id++;
          site_map[std::make_pair(t, bc_t)] = tn;
          outside_sites.push_back(boost::make_tuple(t, cord_t));
        } else {
          tn = site_map.find(std::make_pair(t, bc_t))->second;
        }

        // add bonds
        crossing_bonds.push_back(boost::make_tuple(lattice.index(b), sn, t));
        crossing_bonds.push_back(boost::make_tuple(lattice.index(b), s, tn));
      }
    }
  }

  alps::oxstream oxs;
  oxs << alps::start_tag("GRAPH");

  int id = 1;
  BOOST_FOREACH(const alps::graph_helper<>::site_descriptor& s, lattice.sites()) {
    oxs << alps::start_tag("VERTEX")
        << alps::attribute("id", id++)
        << alps::attribute("type", lattice.site_type(s));
    if (lattice.coordinate(s).begin() != lattice.coordinate(s).end())
      oxs << alps::no_linebreak
          << alps::start_tag("COORDINATE")
          << alps::write_vector(lattice.coordinate(s))
          << alps::end_tag("COORDINATE");
    oxs << alps::end_tag("VERTEX");
  }

  BOOST_FOREACH(site_table_t::value_type const& s, outside_sites) {
    int os = s.get<0>();
    oxs << alps::start_tag("VERTEX")
        << alps::attribute("id", id++)
        << alps::attribute("original_id", os + 1)
        << alps::attribute("type", lattice.site_type(os))
        << alps::attribute("outside", 1);
    if (s.get<1>().size())
      oxs << alps::no_linebreak
          << alps::start_tag("COORDINATE")
          << alps::write_vector(s.get<1>())
          << alps::end_tag("COORDINATE");
    oxs << alps::end_tag("VERTEX");
  }

  id = 1;
  if (unroll) {
    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()) {
      if (!get(alps::boundary_crossing_t(), lattice.graph(), b)) {
        oxs << alps::start_tag("EDGE")
            << alps::attribute("source", lattice.source(b) + 1)
            << alps::attribute("target", lattice.target(b) + 1)
            << alps::attribute("id", id++)
            << alps::attribute("original_id", lattice.index(b) + 1)
            << alps::attribute("type", lattice.bond_type(b));
        if (lattice.bond_vector(b).size())
          oxs << alps::attribute("vector", alps::write_vector(lattice.bond_vector(b)));
        oxs << alps::end_tag("EDGE");
      }
    }
  } else {
    BOOST_FOREACH(const alps::graph_helper<>::bond_descriptor& b, lattice.bonds()) {
      oxs << alps::start_tag("EDGE")
          << alps::attribute("source", lattice.source(b) + 1)
          << alps::attribute("target", lattice.target(b) + 1)
          << alps::attribute("id", lattice.index(b) + 1)
          << alps::attribute("type", lattice.bond_type(b));
      if (lattice.bond_vector(b).size())
        oxs << alps::attribute("vector", alps::write_vector(lattice.bond_vector(b)));
      oxs << alps::end_tag("EDGE");
    }
  }

  BOOST_FOREACH(const bond_table_t::value_type& b, crossing_bonds) {
    alps::graph_helper<>::bond_descriptor ob = lattice.bond(b.get<0>());
    oxs << alps::start_tag("EDGE")
        << alps::attribute("source", b.get<1>() + 1)
        << alps::attribute("target", b.get<2>() + 1)
        << alps::attribute("id", id++)
        << alps::attribute("original_id", lattice.index(ob) + 1)
        << alps::attribute("type", lattice.bond_type(ob))
        << alps::attribute("outside", 1);
    if (lattice.bond_vector(ob).size())
      oxs << alps::attribute("vector", alps::write_vector(lattice.bond_vector(ob)));
    oxs << alps::end_tag("EDGE");
  }

  oxs << alps::end_tag("GRAPH");

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
}
